from django.db import models
import os

class Analysis(models.Model):
    vcf_file = models.FileField(upload_to='uploads/vcf/')
    gff_file = models.FileField(upload_to='uploads/gff/', blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    window_size = models.IntegerField(default=1000, help_text="Tamanho da janela para análise de densidade (bp)")
    
    # Armazenar métricas como JSON para evitar criar muitas colunas
    metrics = models.JSONField(blank=True, null=True)
    
    # Caminhos para imagens geradas (relativo a MEDIA_ROOT)
    plot_quality = models.CharField(max_length=255, blank=True, null=True)
    plot_density = models.CharField(max_length=255, blank=True, null=True)
    
    # Caminho para o CSV de anotação gerado
    annotation_file = models.CharField(max_length=255, blank=True, null=True)

    # Status e Erros
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
