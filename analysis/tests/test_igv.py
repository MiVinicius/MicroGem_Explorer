from django.test import TestCase, Client
from django.urls import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from analysis.models import Analysis
import os

class IGVIntegrationTest(TestCase):
    def setUp(self):
        self.client = Client()
        # Create a dummy VCF file
        self.vcf_content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\t.\t."
        self.vcf_file = SimpleUploadedFile("test.vcf", self.vcf_content, content_type="text/plain")
        
        self.analysis = Analysis.objects.create(vcf_file=self.vcf_file)

    def test_analysis_detail_context_has_vcf_url(self):
        """Test if the analysis detail view provides the VCF URL for IGV."""
        response = self.client.get(reverse('analysis_detail', args=[self.analysis.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Check if the VCF URL is present in the rendered HTML (since we'll put it in JS)
        vcf_url = self.analysis.vcf_file.url
        self.assertContains(response, vcf_url)
        
        # Also check for IGV container div
        self.assertContains(response, 'id="igv-div"')

    def tearDown(self):
        # Clean up files
        if self.analysis.vcf_file:
            if os.path.isfile(self.analysis.vcf_file.path):
                os.remove(self.analysis.vcf_file.path)
