from django.test import TestCase, Client
from django.urls import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from analysis.models import Analysis
import os

class ReferenceGenomeTest(TestCase):
    def setUp(self):
        self.client = Client()
        self.vcf_content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\t.\t."
        self.vcf_file = SimpleUploadedFile("test.vcf", self.vcf_content, content_type="text/plain")
        
        self.fasta_content = b">chr1\nATCG"
        self.fasta_file = SimpleUploadedFile("ref.fasta", self.fasta_content, content_type="text/plain")
        
        self.analysis = Analysis.objects.create(vcf_file=self.vcf_file, reference_file=self.fasta_file)

    def test_igv_uses_custom_reference(self):
        """Test if the IGV configuration uses the custom reference file."""
        response = self.client.get(reverse('analysis_detail', args=[self.analysis.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Check if the reference URL is in the response
        ref_url = self.analysis.reference_file.url
        self.assertContains(response, ref_url)
        
        # Check if 'genome: "hg19"' is NOT present (or handled logic)
        # Actually, if we provide a reference, we might still see hg19 if the logic isn't exclusive, 
        # but we definitely want to see the reference object in the JS options.
        # Let's check for the JS property 'reference:' or 'fastaURL' depending on IGV config
        # IGV.js uses 'reference' object with 'fastaURL'
        
        self.assertContains(response, 'fastaURL')

    def tearDown(self):
        if self.analysis.vcf_file:
            if os.path.isfile(self.analysis.vcf_file.path):
                os.remove(self.analysis.vcf_file.path)
        if self.analysis.reference_file:
            if os.path.isfile(self.analysis.reference_file.path):
                os.remove(self.analysis.reference_file.path)
