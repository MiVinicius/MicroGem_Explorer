from django.core.management.base import BaseCommand
from analysis.vcf_analyzer import VCFAnalyzer
import os

class Command(BaseCommand):
    help = 'Runs QC analysis on a VCF file'

    def add_arguments(self, parser):
        parser.add_argument('--input', type=str, help='Path to the input VCF file', required=True)
        parser.add_argument('--output', type=str, help='Directory to save outputs', default='output')
        parser.add_argument('--window-size', type=int, help='Window size for density analysis', default=1000)
        parser.add_argument('--gff', type=str, help='Path to the GFF file for annotation', required=False)

    def handle(self, *args, **options):
        vcf_path = options['input']
        output_dir = options['output']
        window_size = options['window_size']
        gff_path = options.get('gff')

        self.stdout.write(self.style.SUCCESS(f'Starting QC analysis for {vcf_path}'))

        if not os.path.exists(vcf_path):
            self.stdout.write(self.style.ERROR(f'File not found: {vcf_path}'))
            return

        analyzer = VCFAnalyzer(vcf_path)
        
        try:
            analyzer.parse()
            metrics = analyzer.get_summary()
            
            self.stdout.write(self.style.SUCCESS('QC Metrics:'))
            for key, value in metrics.items():
                self.stdout.write(f'  {key}: {value}')
            
            analyzer.generate_qc_plots(output_dir)
            
            # Density Analysis
            self.stdout.write(f'Calculating mutation density (Window Size: {window_size}bp)...')
            density_df = analyzer.calculate_density(window_size=window_size)
            analyzer.generate_density_plot(density_df, output_dir)
            
            # Annotation
            if gff_path:
                if not os.path.exists(gff_path):
                     self.stdout.write(self.style.ERROR(f'GFF File not found: {gff_path}'))
                else:
                    self.stdout.write(f'Annotating variants using {gff_path}...')
                    annotations = analyzer.annotate_variants(gff_path)
                    if annotations:
                        self.stdout.write(self.style.SUCCESS('Variant Annotations:'))
                        for ann in annotations:
                            self.stdout.write(f"  {ann['CHROM']}:{ann['POS']} ({ann['REF']}->{ann['ALT']}) overlaps: {ann['GENES']}")
                        
                        # Save to CSV
                        import pandas as pd
                        ann_df = pd.DataFrame(annotations)
                        ann_path = os.path.join(output_dir, 'variant_annotations.csv')
                        ann_df.to_csv(ann_path, index=False)
                        self.stdout.write(f"Saved annotations to {ann_path}")
                    else:
                        self.stdout.write("No variants overlapped with genes in the provided GFF.")

            # Generate Summary Report
            analyzer.generate_report(output_dir)

            self.stdout.write(self.style.SUCCESS(f'Analysis complete. Results in {output_dir}'))
            
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Analysis failed: {e}'))
