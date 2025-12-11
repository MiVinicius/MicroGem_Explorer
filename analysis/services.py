import concurrent.futures
import traceback
from .models import Analysis
from .vcf_analyzer import VCFAnalyzer
from django.conf import settings
import os
import pyfaidx

# ThreadPoolExecutor to handle background tasks
executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)

def start_analysis_background(analysis_id):
    """
    Submits the analysis task to the thread pool.
    """
    executor.submit(run_analysis, analysis_id)

def run_analysis(analysis_id):
    """
    Executes the VCF analysis in a separate thread.
    """
    try:
        analysis = Analysis.objects.get(id=analysis_id)
        
        # Update status to PROCESSING
        analysis.status = 'PROCESSING'
        analysis.save()
        
        print(f"Starting background analysis for {analysis_id}...")
        
        # Paths
        vcf_path = analysis.vcf_file.path
        output_rel_dir = f'results/{analysis.id}'
        output_dir = os.path.join(settings.MEDIA_ROOT, output_rel_dir)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Output CSV for table (all variants)
        # We name it 'variants.csv' to distinguish from 'variant_annotations.csv' usually used for GFF
        variants_csv_path = os.path.join(output_dir, 'variants.csv')

        # Run Analysis (Single Pass)
        analyzer = VCFAnalyzer(vcf_path)
        analyzer.process_and_export(variants_csv_path)
        
        # Metrics
        metrics = analyzer.get_summary()
        analysis.metrics = metrics
        
        # Generate Plots (Static backups)
        analyzer.generate_qc_plots(output_dir)
        
        # Density
        # Note: calculate_density now relies on internal state populated by process_and_export
        density_df = analyzer.calculate_density(window_size=analysis.window_size)
        analyzer.generate_density_plot(density_df, output_dir)
        
        # Interactive Data
        plot_data_quality = analyzer.get_quality_distribution_data()
        plot_data_density = analyzer.get_density_data(window_size=analysis.window_size)
        
        analysis.plot_data = {
            'quality': plot_data_quality,
            'density': plot_data_density
        }
        
        # Generator Report
        analyzer.generate_report(output_dir)
        
        # Update paths
        analysis.plot_quality = f'{output_rel_dir}/qc_quality_distribution.png'
        analysis.plot_density = f'{output_rel_dir}/qc_mutation_density.png'
        
        # Annotations (GFF)
        if analysis.gff_file:
            print(f"Annotating variants with GFF: {analysis.gff_file.path}")
            try:
                from .gff_parser import GFFParser
                import pandas as pd
                
                # Parse GFF
                gff_parser = GFFParser(analysis.gff_file.path)
                gff_parser.parse()
                
                # Load variants
                df_variants = pd.read_csv(variants_csv_path)
                
                # Annotate
                annotated_variants = []
                
                # Iterate and annotate
                # Note: This might be slow for many variants/genes. 
                # Optimization: Interval Trees (future work)
                
                # Pre-filter genes by chromosome for speedup
                genes_by_chrom = {}
                for gene in gff_parser.genes:
                    if gene['chrom'] not in genes_by_chrom:
                        genes_by_chrom[gene['chrom']] = []
                    genes_by_chrom[gene['chrom']].append(gene)
                    
                for index, row in df_variants.iterrows():
                    chrom = row['chrom']
                    pos = row['pos']
                    
                    found_genes = []
                    if chrom in genes_by_chrom:
                        for gene in genes_by_chrom[chrom]:
                            if gene['start'] <= pos <= gene['end']:
                                found_genes.append(gene['name'])
                    
                    # Add to list if we want to save separate file or update existing
                    # Let's create a list of dicts for the output CSV
                    var_dict = row.to_dict()
                    var_dict['GENES'] = ";".join(found_genes) if found_genes else "."
                    annotated_variants.append(var_dict)

                # Save Annotated CSV
                annotated_csv_path = os.path.join(output_dir, 'variant_annotations.csv')
                df_annotated = pd.DataFrame(annotated_variants)
                
                # Map columns to uppercase for compatibility with template if needed
                # Template uses: ann.CHROM, ann.POS, ann.REF, ann.ALT, ann.GENES
                # Current df keys are 'chrom', 'pos', etc.
                # Let's rename for the CSV and metrics
                df_annotated.rename(columns={
                    'chrom': 'CHROM', 'pos': 'POS', 'ref': 'REF', 'alt': 'ALT', 
                    'qual': 'QUAL', 'filter': 'FILTER', 'info': 'INFO'
                }, inplace=True)
                
                df_annotated.to_csv(annotated_csv_path, index=False)
                
                # Update Analysis
                analysis.annotation_file = f'{output_rel_dir}/variant_annotations.csv'
                
                # Update metrics with top 50
                # Convert top 50 to dict for JSON serialization
                top_50 = df_annotated.head(50).to_dict('records')
                analysis.metrics['annotations'] = top_50
                
                print(f"Annotation completed. Saved to {annotated_csv_path}")
                
            except Exception as e:
                print(f"Error during annotation: {e}")
                traceback.print_exc()
                # Don't fail the whole analysis just because annotation failed
                analysis.status = 'COMPLETED_WITH_WARNINGS'
        
        # Ensure Reference Index exists
        if analysis.reference_file:
            ref_path = analysis.reference_file.path
            # Check if .fai exists
            fai_path = ref_path + '.fai'
            if not os.path.exists(fai_path):
                print(f"Generating index for {ref_path}...")
                try:
                    # PyFaidx creates index automatically on instanciation if not present
                    pyfaidx.Fasta(ref_path)
                    print(f"Index generated at {fai_path}")
                except Exception as e:
                    print(f"Error generating FASTA index: {e}")
                    # Non-fatal, but IGV might complain

        # Complete
        analysis.status = 'COMPLETED'
        analysis.save()
        print(f"Analysis {analysis_id} completed successfully.")
        
    except Exception as e:
        print(f"Error in background analysis {analysis_id}: {e}")
        traceback.print_exc()
        try:
            analysis = Analysis.objects.get(id=analysis_id)
            analysis.status = 'FAILED'
            analysis.error_message = str(e)
            analysis.save()
        except:
            pass
