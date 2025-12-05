from django.shortcuts import render, redirect, get_object_or_404
from .models import Analysis
from .forms import AnalysisForm
from .vcf_analyzer import VCFAnalyzer
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
            
            # Executar Análise
            vcf_path = analysis.vcf_file.path
            
            # Criar um diretório de saída específico para esta análise
            output_rel_dir = f'results/{analysis.id}'
            output_dir = os.path.join(settings.MEDIA_ROOT, output_rel_dir)
            
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            try:
                analyzer = VCFAnalyzer(vcf_path)
                analyzer.parse()
                
                # Calcular Métricas
                metrics = analyzer.get_summary()
                analysis.metrics = metrics
                
                # Gerar Gráficos (Estáticos - mantidos por compatibilidade ou backup)
                analyzer.generate_qc_plots(output_dir)
                
                density_df = analyzer.calculate_density(window_size=analysis.window_size)
                analyzer.generate_density_plot(density_df, output_dir)
                
                # Pré-calcular dados para gráficos interativos (Plotly) e salvar no banco
                plot_data_quality = analyzer.get_quality_distribution_data()
                plot_data_density = analyzer.get_density_data(window_size=analysis.window_size)
                
                analysis.plot_data = {
                    'quality': plot_data_quality,
                    'density': plot_data_density
                }
                
                # Relatório Resumido
                analyzer.generate_report(output_dir)
                
                # Atualizar modelo com caminhos dos gráficos (relativo a MEDIA_ROOT)
                analysis.plot_quality = f'{output_rel_dir}/qc_quality_distribution.png'
                analysis.plot_density = f'{output_rel_dir}/qc_mutation_density.png'
                
                # Anotação (se GFF fornecido)
                # Anotação (se GFF fornecido)
                if analysis.gff_file:
                    gff_path = analysis.gff_file.path
                    annotations = analyzer.annotate_variants(gff_path)
                    
                    # Salvar para CSV
                    import pandas as pd
                    ann_df = pd.DataFrame(annotations)
                    ann_path = os.path.join(output_dir, 'variant_annotations.csv')
                    ann_df.to_csv(ann_path, index=False)
                    
                    analysis.annotation_file = f'{output_rel_dir}/variant_annotations.csv'

                    # Podemos armazenar anotações em métricas ou em um arquivo/modelo separado
                    # Para simplificar, vamos adicionar os top 50 às métricas para exibição
                    metrics['annotations'] = annotations[:50] 
                
                analysis.save()
                return redirect('analysis_detail', pk=analysis.pk)
                
            except Exception as e:
                # Lidar com erro (excluir modelo ou mostrar mensagem)
                print(f"Erro ao processar análise: {e}")
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
        # Fallback para análises antigas: calcular e salvar
        vcf_path = analysis.vcf_file.path
        try:
            analyzer = VCFAnalyzer(vcf_path)
            analyzer.parse()
            analyzer.calculate_qc_metrics()
            
            plot_data_quality = analyzer.get_quality_distribution_data()
            plot_data_density = analyzer.get_density_data(window_size=analysis.window_size)
            
            # Salvar no banco para a próxima vez
            analysis.plot_data = {
                'quality': plot_data_quality,
                'density': plot_data_density
            }
            analysis.save()
            
        except Exception as e:
            print(f"Erro ao carregar dados de gráfico: {e}")

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

from django.http import JsonResponse

def analysis_variants_api(request, pk):
    """API para retornar variantes paginadas para DataTables."""
    analysis = get_object_or_404(Analysis, pk=pk)
    
    # Parâmetros do DataTables
    start = int(request.GET.get('start', 0))
    length = int(request.GET.get('length', 10))
    draw = int(request.GET.get('draw', 1))
    search_value = request.GET.get('search[value]', None)
    
    vcf_path = analysis.vcf_file.path
    analyzer = VCFAnalyzer(vcf_path)
    
    # Nota: Em produção, o ideal seria não fazer parse() a cada request.
    # Mas como VCFAnalyzer carrega em memória, é o que temos para agora.
    # Uma otimização futura seria cachear ou usar banco de dados.
    
    data, total, filtered = analyzer.get_variants(start, length, search_value)
    
    response = {
        'draw': draw,
        'recordsTotal': total,
        'recordsFiltered': filtered,
        'data': data
    }
    
    return JsonResponse(response)
