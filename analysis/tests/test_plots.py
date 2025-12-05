from django.test import TestCase, Client
from django.urls import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from analysis.models import Analysis
import os
import json

class PlotlyIntegrationTest(TestCase):
    def setUp(self):
        self.client = Client()
        # Create a dummy VCF file
        self.vcf_content = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t30\t.\t.
chr1\t200\t.\tC\tG\t40\t.\t.
chr2\t300\t.\tG\tC\t20\t.\t.
"""
        self.vcf_file = SimpleUploadedFile("test_plots.vcf", self.vcf_content, content_type="text/plain")
        self.analysis = Analysis.objects.create(vcf_file=self.vcf_file)

    def test_analysis_detail_has_plot_data(self):
        """Test if the analysis detail view provides JSON data for Plotly."""
        # We need to trigger the analysis processing first or mock it.
        # Since the view triggers processing on creation, but here we created manually,
        # we might need to call the processing method manually or rely on the view to handle it if it checks for missing metrics.
        # However, the current view logic runs processing ONLY on creation via form.
        # So we should probably manually run the analyzer to populate metrics/plots, 
        # OR update the view to run it if missing.
        # For this test, let's manually populate what we expect the view to pass, 
        # BUT wait, the view calculates context.
        
        # Actually, let's simulate the creation flow to ensure full processing
        # OR just manually run the analyzer logic that the view uses.
        from analysis.vcf_analyzer import VCFAnalyzer
        analyzer = VCFAnalyzer(self.analysis.vcf_file.path)
        analyzer.parse()
        analyzer.calculate_qc_metrics()
        
        # Now that the view is updated, we can test it.
        response = self.client.get(reverse('analysis_detail', args=[self.analysis.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Check if context has plot data
        self.assertIn('plot_data_quality', response.context)
        self.assertIn('plot_data_density', response.context)
        
        # Check if data is valid JSON
        quality_data = json.loads(response.context['plot_data_quality'])
        density_data = json.loads(response.context['plot_data_density'])
        
        self.assertIsInstance(quality_data, list)
        self.assertIsInstance(density_data, dict)

    def tearDown(self):
        if self.analysis.vcf_file:
            if os.path.isfile(self.analysis.vcf_file.path):
                os.remove(self.analysis.vcf_file.path)
