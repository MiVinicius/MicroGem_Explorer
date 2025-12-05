import vcf
import pandas as pd
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import Counter

class VCFAnalyzer:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.reader = None
        self.variants = []
        self.df = None

    def parse(self):
        """Analisa o arquivo VCF e armazena variantes na memória."""
        try:
            self.reader = vcf.Reader(filename=self.vcf_path)
            for record in self.reader:
                self.variants.append(record)
            print(f"Sucesso ao analisar {len(self.variants)} variantes.")
        except Exception as e:
            print(f"Erro ao analisar VCF: {e}")
            raise

    def calculate_qc_metrics(self):
        """Calcula métricas básicas de CQ a partir das variantes analisadas."""
        if not self.variants:
            print("Nenhuma variante para analisar. Chame parse() primeiro.")
            return {}

        qual_scores = [record.QUAL for record in self.variants if record.QUAL is not None]
        
        # Variant Types (SNP vs INDEL)
        # PyVCF record.is_snp is a property
        snp_count = sum(1 for record in self.variants if record.is_snp)
        indel_count = sum(1 for record in self.variants if record.is_indel)
        
        # Filter low quality
        low_qual_count = sum(1 for q in qual_scores if q < 20)
        
        metrics = {
            'total_variants': len(self.variants),
            'snp_count': snp_count,
            'indel_count': indel_count,
            'mean_quality': sum(qual_scores) / len(qual_scores) if qual_scores else 0,
            'low_quality_count': low_qual_count
        }
        
        # Create DataFrame for easier plotting later
        self.df = pd.DataFrame({
            'QUAL': qual_scores
        })
        
        return metrics

    def calculate_ti_tv(self):
        """Calcula a razão Transição/Transversão."""
        if not self.variants:
            return {}

        transitions = 0
        transversions = 0
        
        # Transition: Purine<->Purine (A<->G) or Pyrimidine<->Pyrimidine (C<->T)
        # Transversion: Purine<->Pyrimidine (A<->C, A<->T, G<->C, G<->T)
        
        ti_pairs = {
            ('A', 'G'), ('G', 'A'),
            ('C', 'T'), ('T', 'C')
        }
        
        tv_pairs = {
            ('A', 'C'), ('C', 'A'),
            ('A', 'T'), ('T', 'A'),
            ('G', 'C'), ('C', 'G'),
            ('G', 'T'), ('T', 'G')
        }

        for record in self.variants:
            if not record.is_snp:
                continue
                
            # Handle only the first ALT allele for simplicity in this version
            if len(record.ALT) > 1:
                continue
                
            ref = str(record.REF)
            alt = str(record.ALT[0])
            
            pair = (ref, alt)
            
            if pair in ti_pairs:
                transitions += 1
            elif pair in tv_pairs:
                transversions += 1
                
        ratio = transitions / transversions if transversions > 0 else 0
        
        return {
            'transitions': transitions,
            'transversions': transversions,
            'ti_tv_ratio': ratio
        }

    def generate_qc_plots(self, output_dir):
        """Gera gráficos de CQ e os salva no diretório de saída."""
        if self.df is None or self.df.empty:
            print("Sem dados para plotar.")
            return

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Quality Score Distribution
        # Quality Score Distribution
        plt.figure(figsize=(12, 6))
        # Ensure numeric data
        self.df['QUAL'] = pd.to_numeric(self.df['QUAL'], errors='coerce')
        sns.histplot(self.df['QUAL'].dropna(), kde=True, bins=30)
        plt.title('Distribuição da Pontuação de Qualidade das Variantes')
        plt.xlabel('Pontuação de Qualidade Phred')
        plt.ylabel('Contagem')
        
        output_path = os.path.join(output_dir, 'qc_quality_distribution.png')
        plt.savefig(output_path)
        plt.close()
        print(f"Gráfico de CQ salvo em {output_path}")

    def get_summary(self):
        metrics = self.calculate_qc_metrics()
        ti_tv_metrics = self.calculate_ti_tv()
        metrics.update(ti_tv_metrics)
        return metrics

    def annotate_variants(self, gff_path):
        """Anota variantes com informações de genes de um arquivo GFF."""
        if not self.variants:
            return []

        from analysis.gff_parser import GFFParser
        parser = GFFParser(gff_path)
        parser.parse()

        annotations = []
        for record in self.variants:
            chrom = record.CHROM
            pos = record.POS
            genes = parser.get_genes(chrom, pos)
            
            if genes:
                annotations.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': str(record.REF),
                    'ALT': str(record.ALT[0]) if record.ALT else '.',
                    'GENES': ', '.join(genes)
                })
        
        return annotations

    def calculate_density(self, window_size=1000):
        """Calcula a densidade de variantes em janelas deslizantes."""
        if not self.variants:
            return None

        # Group variants by chromosome
        chrom_variants = {}
        for record in self.variants:
            chrom = record.CHROM
            if chrom not in chrom_variants:
                chrom_variants[chrom] = []
            chrom_variants[chrom].append(record.POS)

        density_data = []

        for chrom, positions in chrom_variants.items():
            if not positions:
                continue
            
            max_pos = max(positions)
            # Create windows
            # Range is inclusive of start, exclusive of end
            for start in range(0, max_pos + window_size, window_size):
                end = start + window_size
                count = sum(1 for pos in positions if start <= pos < end)
                density_data.append({
                    'CHROM': chrom,
                    'WINDOW_START': start,
                    'COUNT': count
                })
        
        return pd.DataFrame(density_data)

    def generate_density_plot(self, density_df, output_dir):
        """Gera gráfico de densidade de mutação."""
        if density_df is None or density_df.empty:
            print("Sem dados de densidade para plotar.")
            return

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        plt.figure(figsize=(12, 6))
        sns.lineplot(data=density_df, x='WINDOW_START', y='COUNT', hue='CHROM', marker='o')
        plt.title('Densidade de Mutação (Janela Deslizante)')
        plt.xlabel('Posição Genômica (bp)')
        plt.ylabel('Contagem de Variantes')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        output_path = os.path.join(output_dir, 'qc_mutation_density.png')
        plt.savefig(output_path)
        plt.close()
        print(f"Gráfico de densidade salvo em {output_path}")

    def generate_report(self, output_dir):
        """Gera um relatório resumido da análise."""
        metrics = self.get_summary()
        
        report_path = os.path.join(output_dir, 'report_summary.txt')
        
        with open(report_path, 'w') as f:
            f.write("MicroGen Variant Explorer - Relatório de Análise\n")
            f.write("================================================\n\n")
            
            f.write("1. Métricas de Controle de Qualidade\n")
            f.write("------------------------------------\n")
            f.write(f"Total de Variantes: {metrics.get('total_variants', 0)}\n")
            f.write(f"Contagem de SNP: {metrics.get('snp_count', 0)}\n")
            f.write(f"Contagem de INDEL: {metrics.get('indel_count', 0)}\n")
            f.write(f"Pontuação Média de Qualidade: {metrics.get('mean_quality', 0):.2f}\n")
            f.write(f"Variantes de Baixa Qualidade (QUAL < 20): {metrics.get('low_quality_count', 0)}\n\n")
            
            f.write("2. Química Genômica (Ti/Tv)\n")
            f.write("---------------------------\n")
            f.write(f"Transições (Ti): {metrics.get('transitions', 0)}\n")
            f.write(f"Transversões (Tv): {metrics.get('transversions', 0)}\n")
            f.write(f"Razão Ti/Tv: {metrics.get('ti_tv_ratio', 0):.2f}\n\n")
            
            f.write("3. Saídas Geradas\n")
            f.write("-----------------\n")
            f.write(f"- Gráfico de Distribuição de Qualidade: {os.path.join(output_dir, 'qc_quality_distribution.png')}\n")
            f.write(f"- Gráfico de Densidade de Mutação: {os.path.join(output_dir, 'qc_mutation_density.png')}\n")
            f.write(f"- Anotações de Variantes: {os.path.join(output_dir, 'variant_annotations.csv')}\n")
            
        print(f"Relatório resumido salvo em {report_path}")
    def get_variants(self, start=0, length=10, search_term=None):
        """
        Retorna uma lista paginada de variantes, opcionalmente filtrada.
        Retorna: (variants_list, total_count, filtered_count)
        """
        # Se variants não estiver carregado, tenta carregar
        if not self.variants:
            try:
                self.parse()
            except:
                return [], 0, 0
        
        # Filtragem
        filtered_variants = self.variants
        if search_term:
            search_term = search_term.lower()
            filtered_variants = [
                v for v in self.variants 
                if search_term in str(v.CHROM).lower() or 
                   search_term in str(v.POS) or
                   search_term in str(v.REF).lower() or
                   (v.ALT and any(search_term in str(a).lower() for a in v.ALT))
            ]
            
        total_count = len(self.variants)
        filtered_count = len(filtered_variants)
        
        # Paginação
        end = start + length
        page_variants = filtered_variants[start:end]
        
        # Formatar para JSON
        data = []
        for v in page_variants:
            data.append({
                'chrom': v.CHROM,
                'pos': v.POS,
                'ref': v.REF,
                'alt': [str(a) for a in v.ALT] if v.ALT else [],
                'qual': v.QUAL,
                'filter': v.FILTER,
                'info': str(v.INFO) # Simplificado para exibição
            })
            
        return data, total_count, filtered_count

    def get_quality_distribution_data(self):
        """Retorna dados da distribuição de qualidade para Plotly."""
        if self.df is None or self.df.empty:
            return []
            
        # Histogram data
        # We can just return the raw values and let Plotly handle the histogram
        # or pre-calculate bins. Returning raw values is more flexible for Plotly.
        return self.df['QUAL'].dropna().tolist()

    def get_density_data(self, window_size=1000):
        """Retorna dados de densidade para Plotly."""
        df = self.calculate_density(window_size)
        if df is None or df.empty:
            return {}
            
        # Format: { 'chrom1': {'x': [pos...], 'y': [count...]}, ... }
        plot_data = {}
        for chrom in df['CHROM'].unique():
            chrom_df = df[df['CHROM'] == chrom]
            plot_data[chrom] = {
                'x': chrom_df['WINDOW_START'].tolist(),
                'y': chrom_df['COUNT'].tolist()
            }
            
        return plot_data
