from django.db import models
import os

class Analysis(models.Model):
    vcf_file = models.FileField(upload_to='uploads/vcf/')
    gff_file = models.FileField(upload_to='uploads/gff/', blank=True, null=True)
    reference_file = models.FileField(
        upload_to='uploads/fasta/',
        blank=True,
        null=True,
        help_text="Arquivo FASTA de referência (opcional)"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    window_size = models.IntegerField(
        default=1000,
        help_text="Tamanho da janela para análise de densidade (bp)"
    )
    
    # Métricas básicas e avançadas armazenadas como JSON
    metrics = models.JSONField(blank=True, null=True)
    
    # Dados pré-calculados para gráficos (qualidade, densidade)
    plot_data = models.JSONField(blank=True, null=True)
    
    # Caminhos para imagens geradas (relativo a MEDIA_ROOT)
    plot_quality = models.CharField(max_length=255, blank=True, null=True)
    plot_density = models.CharField(max_length=255, blank=True, null=True)
    
    # Caminho para CSV de anotações geradas
    annotation_file = models.CharField(max_length=255, blank=True, null=True)

    # Status e mensagens de erro
    STATUS_CHOICES = [
        ('PENDING', 'Pendente'),
        ('PROCESSING', 'Processando'),
        ('COMPLETED', 'Concluído'),
        ('FAILED', 'Falhou'),
    ]
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='PENDING')
    error_message = models.TextField(blank=True, null=True)
    
    def __str__(self):
        return f"Análise {self.id} - {os.path.basename(self.vcf_file.name)}"
    
    # Exemplo de método futuro para contar SNPs ou INDELs
    def count_variant_type(self, variant_type):
        if not self.metrics:
            return 0
        return self.metrics.get(f"{variant_type}_count", 0)
