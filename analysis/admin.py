from django.contrib import admin
from django.utils.html import format_html
from .models import Analysis

@admin.register(Analysis)
class AnalysisAdmin(admin.ModelAdmin):
    list_display = ('id', 'vcf_file', 'created_at', 'status', 'total_variants', 'plot_quality_img')
    list_filter = ('status', 'created_at')
    search_fields = ('vcf_file',)
    readonly_fields = ('metrics', 'plot_quality', 'plot_density', 'annotation_file')

    def total_variants(self, obj):
        if obj.metrics and 'total_variants' in obj.metrics:
            return obj.metrics['total_variants']
        return '-'
    total_variants.short_description = 'Total de Variantes'

    def plot_quality_img(self, obj):
        if obj.plot_quality:
            return format_html('<img src="/media/{}" width="200"/>', obj.plot_quality)
        return '-'
    plot_quality_img.short_description = 'Gráfico QUAL'
