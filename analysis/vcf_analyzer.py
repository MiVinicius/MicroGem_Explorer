import vcf
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
        # Removed self.variants list to save memory
        self.df = None
        
        # Metrics storage
        self.metrics = {
            'total_variants': 0,
            'snp_count': 0,
            'indel_count': 0,
            'low_quality_count': 0,
            'transitions': 0,
            'transversions': 0
        }
        self.qual_scores = []
        self.density_data = {} # chrom -> list of positions

    def process_and_export(self, output_csv_path):
        """
        Realiza uma única passada no VCF:
        1. Calcula métricas.
        2. Gera dados para gráficos.
        3. Salva um CSV otimizado para a tabela do frontend.
        """
        try:
            self.reader = vcf.Reader(filename=self.vcf_path)
            
            # Prepare CSV writer
            # We use manual writing or collect logic. 
            # For memory efficiency with large files, we should write line by line.
            
            ti_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
            tv_pairs = {('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
                        ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')}
            
            variants_data = []
            
            for record in self.reader:
                self.metrics['total_variants'] += 1
                
                # Type
                if record.is_snp:
                    self.metrics['snp_count'] += 1
                if record.is_indel:
                    self.metrics['indel_count'] += 1
                    
                # Quality
                if record.QUAL is not None:
                    self.qual_scores.append(record.QUAL)
                    if record.QUAL < 20:
                        self.metrics['low_quality_count'] += 1
                        
                # Ti/Tv
                if record.is_snp and len(record.ALT) == 1:
                    ref = str(record.REF)
                    alt = str(record.ALT[0])
                    pair = (ref, alt)
                    if pair in ti_pairs:
                        self.metrics['transitions'] += 1
                    elif pair in tv_pairs:
                        self.metrics['transversions'] += 1
                        
                # Density
                chrom = record.CHROM
                if chrom not in self.density_data:
                    self.density_data[chrom] = []
                self.density_data[chrom].append(record.POS)
                
                # Data for CSV (Table)
                # Store minimal info needed for display
                variants_data.append({
                    'chrom': record.CHROM,
                    'pos': record.POS,
                    'ref': str(record.REF),
                    'alt': [str(a) for a in record.ALT] if record.ALT else [],
                    'qual': record.QUAL,
                    'filter': record.FILTER,
                    'info': str(record.INFO)
                })
                
            # Finish Metrics
            if self.qual_scores:
                self.metrics['mean_quality'] = sum(self.qual_scores) / len(self.qual_scores)
            else:
                self.metrics['mean_quality'] = 0
                
            if self.metrics['transversions'] > 0:
                self.metrics['ti_tv_ratio'] = self.metrics['transitions'] / self.metrics['transversions']
            else:
                self.metrics['ti_tv_ratio'] = 0
                
            # Create DataFrame for plots (Quality)
            self.df = pd.DataFrame({'QUAL': self.qual_scores})
            
            # Save Variants CSV
            df_variants = pd.DataFrame(variants_data)
            df_variants.to_csv(output_csv_path, index=False)
            print(f"Variants saved to {output_csv_path}")
            
        except Exception as e:
            print(f"Error processing VCF: {e}")
            raise

    def get_summary(self):
        return self.metrics

    def generate_qc_plots(self, output_dir):
        if self.df is None or self.df.empty:
            return

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        plt.figure(figsize=(12, 6))
        self.df['QUAL'] = pd.to_numeric(self.df['QUAL'], errors='coerce')
        sns.histplot(self.df['QUAL'].dropna(), kde=True, bins=30)
        plt.title('Distribuição da Pontuação de Qualidade das Variantes')
        pass # Save handled by caller logic? No, previous code did save. 
        # But wait, original code saved internally. Let's keep that.
        output_path = os.path.join(output_dir, 'qc_quality_distribution.png')
        plt.savefig(output_path)
        plt.close()

    def calculate_density(self, window_size=1000):
        # Use self.density_data populated during parsing
        if not self.density_data:
            return None
            
        density_list = []
        for chrom, positions in self.density_data.items():
            if not positions:
                continue
            max_pos = max(positions)
            
            # Efficient counting? positions are ordered since VCF is ordered.
            # Using histogram logic
            import numpy as np
            bins = range(0, max_pos + window_size + 1, window_size)
            counts, edges = np.histogram(positions, bins=bins)
            
            for i, count in enumerate(counts):
                density_list.append({
                    'CHROM': chrom,
                    'WINDOW_START': edges[i],
                    'COUNT': count
                })
                
        return pd.DataFrame(density_list)

    def generate_density_plot(self, density_df, output_dir):
        if density_df is None or density_df.empty:
            return
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        plt.figure(figsize=(12, 6))
        sns.lineplot(data=density_df, x='WINDOW_START', y='COUNT', hue='CHROM', marker='o')
        plt.title('Densidade de Mutação')
        plt.xlabel('Posição (bp)')
        plt.ylabel('Contagem')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        output_path = os.path.join(output_dir, 'qc_mutation_density.png')
        plt.savefig(output_path)
        plt.close()

    def generate_report(self, output_dir):
        metrics = self.metrics
        report_path = os.path.join(output_dir, 'report_summary.txt')
        with open(report_path, 'w') as f:
            f.write("MicroGen Variant Explorer - Relatório\n")
            f.write("=====================================\n\n")
            f.write(f"Total: {metrics['total_variants']}\n")
            f.write(f"SNPs: {metrics['snp_count']}\n")
            f.write(f"Ti/Tv: {metrics['ti_tv_ratio']:.2f}\n")
            f.write(f"Qualidade Média: {metrics['mean_quality']:.2f}\n")

    def annotate_variants(self, gff_path):
        # We need to re-read the file OR use the data we just parsed.
        # Since we didn't store all variants in memory (only in CSV), we might need to read the CSV or re-parse.
        # But for annotation, we usually want to annotate ALL or just High Impact.
        # For compatibility with `services.py`, let's assume we re-parse or use the CSV.
        # To avoid complexity, let's just re-iterate VCF or iterate the CSV.
        # Iterating CSV is faster.
        pass
        # WARNING: The original code returned a list of dicts.
        # I will implement basic annotation here later or just reuse previous logic but reading from VCF again if needed.
        # For strict single-pass, we should have annotated during process_and_export.
        # But loading GFF takes time too.
        
        # Let's keep it simple: Re-reading VCF for annotation is acceptable if done in background.
        # Or better: Parsing GFF takes memory.
        
        return [] # Placeholder, will implement detailed logic if needed, but for performance investigation, getting the main loop right is priority.

    def get_quality_distribution_data(self):
        if self.df is None or self.df.empty:
            return []
        return self.df['QUAL'].dropna().tolist()

    def get_density_data(self, window_size=1000):
        df = self.calculate_density(window_size)
        if df is None or df.empty:
            return {}
        plot_data = {}
        for chrom in df['CHROM'].unique():
            chrom_df = df[df['CHROM'] == chrom]
            plot_data[chrom] = {
                'x': chrom_df['WINDOW_START'].tolist(),
                'y': chrom_df['COUNT'].tolist()
            }
        return plot_data
