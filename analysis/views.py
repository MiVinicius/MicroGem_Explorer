from django.shortcuts import render, redirect, get_object_or_404
from .models import Analysis
from .forms import AnalysisForm
from .vcf_analyzer import VCFAnalyzer
from .services import start_analysis_background
import os
from django.conf import settings

def home(request):
    return render(request, 'analysis/home.html')

def analysis_list(request):
    analyses = Analysis.objects.all().order_by('-created_at')
    return render(request, 'analysis/analysis_list.html', {'analyses': analyses})

def analysis_create(request):
    if request.method == 'POST':
        form = AnalysisForm(request.POST, request.FILES)
        if form.is_valid():
            analysis = form.save()
            
            try:
                # Iniciar análise em background
                start_analysis_background(analysis.id)
                return redirect('analysis_detail', pk=analysis.pk)
                
            except Exception as e:
                # Lidar com erro (excluir modelo ou mostrar mensagem)
                print(f"Erro ao iniciar análise: {e}")
                analysis.delete()
                form.add_error(None, f"Erro ao processar arquivo: {e}")
    else:
        form = AnalysisForm()
    return render(request, 'analysis/analysis_form.html', {'form': form})

def analysis_detail(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    
    # Arredondar métricas para exibição para evitar problemas com filtros de template
    if analysis.metrics:
        if 'mean_quality' in analysis.metrics:
            analysis.metrics['mean_quality'] = round(analysis.metrics['mean_quality'], 2)
        if 'ti_tv_ratio' in analysis.metrics:
            analysis.metrics['ti_tv_ratio'] = round(analysis.metrics['ti_tv_ratio'], 2)
            
    # Obter dados para gráficos interativos
    import json
    
    plot_data_quality = []
    plot_data_density = {}
    
    # Tentar usar dados em cache
    if analysis.plot_data:
        plot_data_quality = analysis.plot_data.get('quality', [])
        plot_data_density = analysis.plot_data.get('density', {})
    else:
        # Se não houver dados, não tentamos calcular síncrono para evitar travar.
        # O background service irá preencher isso.
        pass

    context = {
        'analysis': analysis,
        'plot_data_quality': json.dumps(plot_data_quality),
        'plot_data_density': json.dumps(plot_data_density)
    }
            
    return render(request, 'analysis/analysis_detail.html', context)

def analysis_delete(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    if request.method == 'POST':
        analysis.delete()
        return redirect('analysis_list')
    return render(request, 'analysis/analysis_confirm_delete.html', {'analysis': analysis})
    return render(request, 'analysis/analysis_confirm_delete.html', {'analysis': analysis})

import pandas as pd
from django.http import JsonResponse

def analysis_variants_api(request, pk):
    """API para retornar variantes paginadas para DataTables."""
    analysis = get_object_or_404(Analysis, pk=pk)
    
    # Parâmetros do DataTables
    start = int(request.GET.get('start', 0))
    length = int(request.GET.get('length', 10))
    draw = int(request.GET.get('draw', 1))
    search_value = request.GET.get('search[value]', None)
    
    # Caminho do CSV gerado pelo services.py
    output_rel_dir = f'results/{analysis.id}'
    csv_path = os.path.join(settings.MEDIA_ROOT, output_rel_dir, 'variants.csv')
    
    data = []
    total = 0
    filtered = 0
    
    if os.path.exists(csv_path):
        try:
            # Ler CSV otimizado
            # Para arquivos muito grandes, usar chunks. Aqui vamos ler tudo pois é mais rápido que parsear VCF
            # Se for > 100MB, devemos usar chunks ou banco de dados.
            df = pd.read_csv(csv_path)
            
            total = len(df)
            
            # Busca
            if search_value:
                df = df[
                    df.astype(str).apply(lambda x: x.str.contains(search_value, case=False, na=False)).any(axis=1)
                ]
            
            filtered = len(df)
            
            # Paginação
            df_page = df.iloc[start:start+length]
            
            # Converter para dict
            data = df_page.to_dict('records')
            
        except Exception as e:
            print(f"Erro ao ler CSV de variantes: {e}")
    
    response = {
        'draw': draw,
        'recordsTotal': total,
        'recordsFiltered': filtered,
        'data': data
    }
    
    return JsonResponse(response)
