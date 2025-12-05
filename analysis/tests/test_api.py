from django.test import TestCase, Client
from django.urls import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from analysis.models import Analysis
import json
import os

class VariantAPITest(TestCase):
    def setUp(self):
        self.client = Client()
        # Create a dummy VCF file with a few variants
        self.vcf_content = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t30\t.\t.
chr1\t200\t.\tC\tG\t40\t.\t.
chr2\t300\t.\tG\tC\t50\t.\t.
"""
        self.vcf_file = SimpleUploadedFile("test_api.vcf", self.vcf_content, content_type="text/plain")
        self.analysis = Analysis.objects.create(vcf_file=self.vcf_file)

    def test_api_pagination(self):
        """Test if the API returns paginated results."""
        url = reverse('analysis_variants_api', args=[self.analysis.pk])
        
        # Request first page, size 1
        response = self.client.get(url, {'start': 0, 'length': 1, 'draw': 1})
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        self.assertEqual(data['draw'], 1)
        self.assertEqual(data['recordsTotal'], 3)
        self.assertEqual(data['recordsFiltered'], 3)
        self.assertEqual(len(data['data']), 1)
        self.assertEqual(data['data'][0]['pos'], 100)

        # Request second page
        response = self.client.get(url, {'start': 1, 'length': 1, 'draw': 2})
        data = json.loads(response.content)
        self.assertEqual(len(data['data']), 1)
        self.assertEqual(data['data'][0]['pos'], 200)

    def test_api_search(self):
        """Test if the API filters results by search term (CHROM)."""
        url = reverse('analysis_variants_api', args=[self.analysis.pk])
        
        # Search for "chr2"
        response = self.client.get(url, {'start': 0, 'length': 10, 'draw': 1, 'search[value]': 'chr2'})
        data = json.loads(response.content)
        
        self.assertEqual(data['recordsTotal'], 3)
        self.assertEqual(data['recordsFiltered'], 1)
        self.assertEqual(len(data['data']), 1)
        self.assertEqual(data['data'][0]['chrom'], 'chr2')

    def tearDown(self):
        if self.analysis.vcf_file:
            if os.path.isfile(self.analysis.vcf_file.path):
                os.remove(self.analysis.vcf_file.path)
