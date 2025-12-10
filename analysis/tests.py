from django.test import TestCase
from django.core.files.uploadedfile import SimpleUploadedFile
from .models import Analysis
from .services import run_analysis
import os

class AnalysisModelTest(TestCase):
    def setUp(self):
        # Cria um arquivo VCF simples para teste
        self.vcf_file = SimpleUploadedFile(
            "test.vcf",
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tT\t50\tPASS\t.\n",
            content_type="text/plain"
        )
        
        self.analysis = Analysis.objects.create(
            vcf_file=self.vcf_file,
            window_size=1000
        )

    def test_analysis_creation(self):
        self.assertEqual(self.analysis.status, 'PENDING')
        self.assertIsNone(self.analysis.metrics)

class AnalysisServiceTest(TestCase):
    def setUp(self):
        self.vcf_file = SimpleUploadedFile(
            "test.vcf",
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tT\t50\tPASS\t.\n",
            content_type="text/plain"
        )
        self.analysis = Analysis.objects.create(
            vcf_file=self.vcf_file,
            window_size=1000
        )

    def test_run_analysis_background(self):
        # Executa a análise de forma síncrona para teste
        run_analysis(self.analysis.id)
        analysis = Analysis.objects.get(id=self.analysis.id)
        self.assertEqual(analysis.status, 'COMPLETED')
        self.assertIn('total_variants', analysis.metrics)
        self.assertEqual(analysis.metrics['total_variants'], 1)
