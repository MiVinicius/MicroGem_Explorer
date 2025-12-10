import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from collections import defaultdict, Counter
import vcf  # PyVCF
from BCBio import GFF
import warnings
warnings.filterwarnings("ignore")


class VCFAnalyzer:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.df_variants = None
        self.df_quality = None
        self.qual_scores = []
        self.density_data = defaultdict(list)
        self.annotations = []

        self.metrics = {
            "total_variants": 0,
            "snp_count": 0,
            "indel_count": 0,
            "mnv_count": 0,
            "low_quality_count": 0,
            "transitions": 0,
            "transversions": 0,
            "ti_tv_ratio": 0,
            "mean_quality": 0,
            "chrom_distribution": {},
            "annotations": [],
            "impact_counts": defaultdict(int),
            "top_genes": [],
            "ti_tv_gene": {},
            "hotspots": []
        }

    # ---------------------------
    # PROCESSAMENTO DO VCF
    # ---------------------------
    def process_and_export(self, output_csv_path):
        reader = vcf.Reader(filename=self.vcf_path)
        variants_list = []

        TI = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        TV = {("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
              ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}

        for record in reader:
            self.metrics["total_variants"] += 1
            chrom = record.CHROM
            pos = record.POS
            ref = str(record.REF)
            alt_list = [str(a) for a in record.ALT] if record.ALT else ["N"]
            qual = record.QUAL if record.QUAL is not None else 0

            var_type = "SNP" if record.is_snp else "INDEL" if record.is_indel else "MNV"
            if record.is_snp:
                self.metrics["snp_count"] += 1
            elif record.is_indel:
                self.metrics["indel_count"] += 1
            else:
                self.metrics["mnv_count"] += 1

            self.qual_scores.append(qual)
            if qual < 20:
                self.metrics["low_quality_count"] += 1

            if record.is_snp and len(alt_list) == 1:
                pair = (ref, alt_list[0])
                if pair in TI:
                    self.metrics["transitions"] += 1
                elif pair in TV:
                    self.metrics["transversions"] += 1

            self.metrics["chrom_distribution"][chrom] = self.metrics["chrom_distribution"].get(chrom, 0) + 1
            self.density_data[chrom].append(pos if pos else 0)

            # Impacto funcional
            impact_info = record.INFO.get("ANN")
            if impact_info:
                ann_str = impact_info[0] if isinstance(impact_info, list) else impact_info
                if "synonymous_variant" in ann_str:
                    self.metrics["impact_counts"]["synonymous"] += 1
                elif "missense_variant" in ann_str:
                    self.metrics["impact_counts"]["nonsynonymous"] += 1
                elif "frameshift_variant" in ann_str:
                    self.metrics["impact_counts"]["frameshift"] += 1
                elif "stop_gained" in ann_str:
                    self.metrics["impact_counts"]["stop_gain"] += 1
                elif "stop_lost" in ann_str:
                    self.metrics["impact_counts"]["stop_loss"] += 1

            variants_list.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": ",".join(alt_list),
                "QUAL": qual,
                "TYPE": var_type,
                "GENES": ""
            })

        self.qual_scores = [q for q in self.qual_scores if q is not None]
        if self.qual_scores:
            self.metrics["mean_quality"] = float(np.mean(self.qual_scores))
        if self.metrics["transversions"] > 0:
            self.metrics["ti_tv_ratio"] = self.metrics["transitions"] / self.metrics["transversions"]

        df = pd.DataFrame(variants_list)
        df.to_csv(output_csv_path, index=False)
        self.df_variants = df
        self.df_quality = pd.DataFrame({"QUAL": self.qual_scores})

        if not self.qual_scores:
            self.qual_scores = [0]
        if not self.density_data:
            self.density_data = {"chr1": [0]}

        return df

    # ---------------------------
    # ANOTAÇÃO COM GFF
    # ---------------------------
    def annotate_with_gff(self, gff_path):
        if self.df_variants is None:
            raise ValueError("VCF não processado")

        annotations = []
        gff_data = {}
        with open(gff_path) as gff_file:
            for rec in GFF.parse(gff_file):
                chrom = rec.id
                gff_data.setdefault(chrom, [])
                for feature in rec.features:
                    if 'gene' in feature.type.lower():
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        gene_id = feature.qualifiers.get('Name', feature.qualifiers.get('ID', ['unknown']))[0]
                        gff_data[chrom].append((start, end, gene_id))

        for idx, row in self.df_variants.iterrows():
            chrom = row['CHROM']
            pos = row['POS']
            genes = []
            if chrom in gff_data:
                for start, end, gene_id in gff_data[chrom]:
                    if start <= pos <= end:
                        genes.append(gene_id)
            self.df_variants.at[idx, 'GENES'] = ','.join(genes) if genes else ''
            annotations.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": row['REF'],
                "ALT": row['ALT'],
                "GENES": ','.join(genes) if genes else ''
            })

        self.annotations = annotations
        self.metrics['annotations'] = annotations

        # Top genes
        gene_counter = Counter()
        for var in annotations:
            for gene in var.get("GENES", "").split(","):
                if gene:
                    gene_counter[gene] += 1
        self.metrics["top_genes"] = gene_counter.most_common(10)

        # Ti/Tv por gene
        TI = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        TV = {("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
              ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}
        ti_tv_gene = defaultdict(lambda: {"Ti": 0, "Tv": 0})
        for var in annotations:
            genes = [g for g in var.get("GENES", "").split(",") if g]
            ref, alt = var["REF"], var["ALT"]
            for g in genes:
                if (ref, alt) in TI:
                    ti_tv_gene[g]["Ti"] += 1
                elif (ref, alt) in TV:
                    ti_tv_gene[g]["Tv"] += 1
        self.metrics["ti_tv_gene"] = {g: t["Ti"] / t["Tv"] if t["Tv"] > 0 else 0
                                      for g, t in ti_tv_gene.items()}

        return annotations

    # ---------------------------
    # DENSIDADE & HOTSPOTS
    # ---------------------------
    def calculate_density(self, window_size=1000):
        density_rows = []
        for chrom, positions in self.density_data.items():
            if not positions:
                continue
            max_pos = max(positions)
            bins = range(0, max_pos + window_size, window_size)
            counts, edges = np.histogram(positions, bins=bins)
            for i, count in enumerate(counts):
                density_rows.append({
                    "CHROM": chrom,
                    "WINDOW_START": int(edges[i]),
                    "COUNT": int(count),
                    "DENSITY_NORM": float(count / (window_size / 1000))  # var/kb
                })
        df = pd.DataFrame(density_rows)
        self.metrics["hotspots"] = df[df["DENSITY_NORM"] > 5].to_dict(orient="records")
        return df

    def get_density_data(self, window_size=1000):
        df = self.calculate_density(window_size)
        if df.empty:
            return {}
        plot_data = {}
        for chrom in df["CHROM"].unique():
            sub = df[df["CHROM"] == chrom]
            plot_data[chrom] = {
                "x": sub["WINDOW_START"].tolist(),
                "y": sub["DENSITY_NORM"].tolist(),
                "count": sub["COUNT"].tolist()
            }
        return plot_data

    # ---------------------------
    # PLOTS
    # ---------------------------
    def generate_qc_plots(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        # Gráfico QUAL
        qual_path = os.path.join(output_dir, "qc_quality_distribution.png")
        if self.df_quality is not None:
            plt.figure(figsize=(8, 4))
            plt.hist(self.df_quality["QUAL"], bins=30, color='green', alpha=0.7)
            plt.xlabel("QUAL")
            plt.ylabel("Contagem")
            plt.title("Distribuição de Qualidade")
            plt.tight_layout()
            plt.savefig(qual_path)
            plt.close()

        # Gráfico Densidade por cromossomo
        density_path = os.path.join(output_dir, "density_per_chrom.png")
        density_data = self.get_density_data()
        if density_data:
            plt.figure(figsize=(10, 5))
            for chrom, data in density_data.items():
                plt.bar(data["x"], data["y"], width=900, alpha=0.6, label=chrom)
            plt.xlabel("Posição Genômica (bp)")
            plt.ylabel("Densidade de Variantes (var/kb)")
            plt.title("Densidade de Variantes por Cromossomo")
            plt.legend()
            plt.tight_layout()
            plt.savefig(density_path)
            plt.close()

        return {"qual_plot": qual_path, "density_plot": density_path}

    # ---------------------------
    # MÉTRICAS
    # ---------------------------
    def get_summary(self):
        return self.metrics

    def get_quality_distribution_data(self):
        return self.qual_scores
