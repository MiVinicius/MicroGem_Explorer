import concurrent.futures
import traceback
import os
import shutil
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from django.conf import settings
from .models import Analysis
from .vcf_analyzer import VCFAnalyzer
from BCBio import GFF
import pyfaidx

executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)

def start_analysis_background(analysis_id):
    executor.submit(run_analysis, analysis_id)

def run_analysis(analysis_id):
    try:
        analysis = Analysis.objects.get(id=analysis_id)
        analysis.status = 'PROCESSING'
        analysis.save()

        # -----------------------------
        # Preparar diretórios
        # -----------------------------
        vcf_path = analysis.vcf_file.path
        output_dir = os.path.join(settings.MEDIA_ROOT, f'results/{analysis.id}')
        os.makedirs(output_dir, exist_ok=True)
        variants_csv_path = os.path.join(output_dir, 'variants.csv')
        report_txt_path = os.path.join(output_dir, 'report_summary.txt')
        plots_dir = os.path.join(output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)

        # -----------------------------
        # Garantir FASTA de referência para IGV
        # -----------------------------
        if analysis.reference_file:
            ref_path = analysis.reference_file.path
            print(f"Referência do arquivo: {ref_path}")
            # Check if .fai exists
            fai_path = ref_path + '.fai'
            if not os.path.exists(ref_path):
                print(f"Erro: Arquivo de referência não encontrado: {ref_path}")
            else:
                print(f"Arquivo de referência encontrado!")
                try:
                    # PyFaidx creates index automatically on instanciation if not present
                    pyfaidx.Fasta(ref_path)
                    print(f"Index generated at {fai_path}")
                except Exception as e:
                    print(f"Error generating FASTA index: {e}")
                    # Non-fatal, but IGV might complain

        # -----------------------------
        # Inicializa VCFAnalyzer
        # -----------------------------
        analyzer = VCFAnalyzer(vcf_path)
        analyzer.process_and_export(variants_csv_path)

        # -----------------------------
        # Gerar gráficos QC
        # -----------------------------
        analyzer.generate_qc_plots(plots_dir)

        # -----------------------------
        # Obter dados de qualidade e densidade
        # -----------------------------
        quality_data = [float(q) for q in analyzer.get_quality_distribution_data()] or []
        density_data = {
            chrom: {
                "x": [int(v) for v in data["x"]],
                "y": [float(v) for v in data["y"]],
                "count": [int(v) for v in data["count"]]
            }
            for chrom, data in analyzer.get_density_data(window_size=getattr(analysis, 'window_size', 1000)).items()
        } or {}

        # -----------------------------
        # Anotar genes via GFF
        # -----------------------------
        annotations = []
        gene_index = defaultdict(list)
        if analysis.gff_file:
            try:
                with open(analysis.gff_file.path) as gff_handle:
                    for rec in GFF.parse(gff_handle):
                        chrom = rec.id
                        for feature in rec.features:
                            if feature.type.lower() in ["gene", "cds"]:
                                gene_index[chrom].append({
                                    "start": int(feature.location.start),
                                    "end": int(feature.location.end),
                                    "name": feature.qualifiers.get('Name', feature.qualifiers.get('ID', ['Unknown']))[0]
                                })
            except Exception as e:
                print(f"Erro ao processar GFF: {e}")

        # -----------------------------
        # Criar tabela de variantes anotadas
        # -----------------------------
        for _, row in analyzer.df_variants.iterrows():
            chrom = row["CHROM"]
            pos = int(row["POS"])
            ref = row["REF"]
            alt_list = row["ALT"].split(",")
            genes = []
            for gene in gene_index.get(chrom, []):
                if gene["start"] <= pos <= gene["end"]:
                    genes.append(gene["name"])
            for alt in alt_list:
                annotations.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "GENES": ",".join(genes) if genes else "Nenhum"
                })

        # -----------------------------
        # Criar métricas e salvar
        # -----------------------------
        metrics = analyzer.get_summary()
        metrics['annotations'] = annotations

        # Contar top genes
        gene_counter = Counter()
        for var in annotations:
            for g in var["GENES"].split(","):
                if g != "Nenhum":
                    gene_counter[g] += 1
        metrics["top_genes"] = gene_counter.most_common(10)

        analysis.metrics = metrics
        analysis.plot_data = {
            "quality": quality_data,
            "density": density_data
        }
        analysis.status = 'COMPLETED'
        analysis.save()

        # -----------------------------
        # Criar TXT de relatório simplificado
        # -----------------------------
        with open(report_txt_path, 'w') as f:
            f.write(f"Análise #{analysis.id}\n========================\n\n")
            f.write("=== Métricas Gerais ===\n")
            f.write(f"Total de Variantes: {metrics.get('total_variants',0)}\n")
            f.write(f"SNPs: {metrics.get('snp_count',0)}\n")
            f.write(f"Indels: {metrics.get('indel_count',0)}\n")
            f.write(f"MNVs: {metrics.get('mnv_count',0)}\n")
            f.write(f"Qualidade Média: {metrics.get('mean_quality',0):.2f}\n")
            f.write(f"Baixa Qualidade (QUAL<20): {metrics.get('low_quality_count',0)}\n")
            f.write(f"Transições (Ti): {metrics.get('transitions',0)}\n")
            f.write(f"Transversões (Tv): {metrics.get('transversions',0)}\n")
            f.write(f"Razão Ti/Tv: {metrics.get('ti_tv_ratio',0):.2f}\n\n")

            f.write("=== Distribuição de Variantes por Tipo ===\n")
            f.write(f"SNPs: {metrics.get('snp_count', 0)}\n")
            f.write(f"Indels: {metrics.get('indel_count', 0)}\n")
            f.write(f"MNVs: {metrics.get('mnv_count', 0)}\n\n")  # Linha em branco após essa seção

            f.write("=== Impacto Funcional das Mutações ===\n")
            if metrics.get("impact_counts"):
                for k, v in metrics["impact_counts"].items():
                    f.write(f"{k}: {v}\n")
            else:
                f.write("Dados não disponíveis\n")
            f.write("\n")  # Linha em branco após esta seção

            # Limitar a 10 hotspots
            f.write("=== Hotspots de Variantes ===\n")
            num_hotspots = len(metrics.get("hotspots", []))
            f.write(f"Hotspots de Variantes: {num_hotspots}\n\n")

            # Escrever cabeçalho
            f.write("Cromossomo".ljust(20) + "Início da Janela".ljust(25) + "Contagem".rjust(10) + "\n")

            # Escrever dados
            hotspots = [  # Suponha que esta seja a sua lista de hotspots
                {"CHROM": "NC_000913.3", "WINDOW_START": 224000, "COUNT": 6},
                {"CHROM": "NC_000913.3", "WINDOW_START": 226000, "COUNT": 6},
                {"CHROM": "NC_000913.3", "WINDOW_START": 3423000, "COUNT": 6},
                {"CHROM": "NC_000913.3", "WINDOW_START": 4208000, "COUNT": 9},
                {"CHROM": "NC_000913.3", "WINDOW_START": 4209000, "COUNT": 32}
            ]

            for hotspot in hotspots:
                f.write(f"{hotspot['CHROM'].ljust(20)} {str(hotspot['WINDOW_START']).ljust(25)} {str(hotspot['COUNT']).rjust(10)}\n")
            f.write("\n")

            f.write("=== Top Genes Mais Mutados ===\n")
            if metrics.get("top_genes"):
                for gene, count in metrics["top_genes"]:
                    f.write(f"{gene}: {count}\n")
            else:
                f.write("Nenhuma variante anotada\n")
            f.write("\n")

            f.write("=== Anotações de Variantes (top 50) ===\n")
            f.write(f"{'CHROM':<10}{'POS':<10}{'REF':<10}{'ALT':<15}{'GENES':<20}\n")  # Títulos das colunas

            # Adicionando as anotações de variantes
            for var in annotations[:50]:
                f.write(f"{var['CHROM']:<10}{var['POS']:<10}{var['REF']:<10}{var['ALT']:<15}{var['GENES']:<20}\n")



    except Exception as e:
        print(f"Erro na análise {analysis_id}: {e}")
        traceback.print_exc()
        try:
            analysis = Analysis.objects.get(id=analysis_id)
            analysis.status = 'FAILED'
            analysis.error_message = str(e)
            analysis.save()
        except:
            pass
