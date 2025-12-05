from django import forms
from .models import Analysis

class AnalysisForm(forms.ModelForm):
    class Meta:
        model = Analysis
        fields = ['vcf_file', 'gff_file', 'reference_file', 'window_size']
        widgets = {
            'vcf_file': forms.FileInput(attrs={'class': 'form-control'}),
            'gff_file': forms.FileInput(attrs={'class': 'form-control'}),
            'reference_file': forms.FileInput(attrs={'class': 'form-control'}),
            'window_size': forms.NumberInput(attrs={'class': 'form-control'}),
        }
