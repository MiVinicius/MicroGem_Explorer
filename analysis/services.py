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
            # For now, we skip or simple implementation. 
            # In optimized version, we might re-read the generated CSV to annotate?
            # Or re-read VCF. Let's look at skipping to ensure basic flow works first.
            # If user needs annotation, we can add it later.
            pass
        
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
