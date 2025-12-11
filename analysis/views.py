from django.shortcuts import render, redirect, get_object_or_404
from .models import Analysis
from .forms import AnalysisForm
from .services import start_analysis_background
from django.conf import settings
import os
import pandas as pd
import json
from django.http import JsonResponse

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
                # Inicia análise em background
                start_analysis_background(analysis.id)
                return redirect('analysis_detail', pk=analysis.pk)
            except Exception as e:
                analysis.delete()
                form.add_error(None, f"Erro ao processar arquivo: {e}")
    else:
        form = AnalysisForm()
    return render(request, 'analysis/analysis_form.html', {'form': form})

def analysis_detail(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)

    # Métricas
    metrics = analysis.metrics or {}

    # --- ARREDONDAR MÉTRICAS ---
    for key in ['mean_quality', 'ti_tv_ratio']:
        if key in metrics and metrics[key] is not None:
            try:
                metrics[key] = round(float(metrics[key]), 2)
            except:
                pass
    analysis.metrics = metrics



    # ---- Separar dados ----

    # Variantes (vem do CSV)
    variants = []
    variants_csv = os.path.join(settings.MEDIA_ROOT, f"results/{analysis.id}/variants.csv")

    if os.path.exists(variants_csv):
        try:
            df = pd.read_csv(variants_csv)
            variants = df.fillna("").astype(str).to_dict("records")
        except Exception as e:
            print("Erro ao carregar variants.csv:", e)

    # Anotações (vem das métricas JSON)
    annotations = metrics.get("annotations", [])

    context = {
        'analysis': analysis,
        'analysis_metrics': metrics,
        'variants_json': json.dumps(variants),
        'annotations_json': json.dumps(annotations),
        'MEDIA_URL': settings.MEDIA_URL,
    }

    return render(request, 'analysis/analysis_detail.html', context)


def analysis_delete(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    if request.method == 'POST':
        analysis.delete()
        return redirect('analysis_list')
    return render(request, 'analysis/analysis_confirm_delete.html', {'analysis': analysis})

def analysis_variants_api(request, pk):
    """API para DataTables retornar variantes paginadas."""
    analysis = get_object_or_404(Analysis, pk=pk)

    start = int(request.GET.get('start', 0))
    length = int(request.GET.get('length', 10))
    draw = int(request.GET.get('draw', 1))
    search_value = request.GET.get('search[value]', '')

    csv_path = os.path.join(settings.MEDIA_ROOT, f'results/{analysis.id}/variants.csv')

    data = []
    records_total = 0
    records_filtered = 0

    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path, encoding='utf-8')
            records_total = len(df)

            if search_value:
                df = df[df.astype(str).apply(lambda x: x.str.contains(search_value, case=False, na=False)).any(axis=1)]

            records_filtered = len(df)
            df_page = df.iloc[start:start+length]
            data = df_page.astype(str).to_dict('records')

        except Exception as e:
            print(f"Erro ao ler CSV de variantes: {e}")

    return JsonResponse({
        'draw': draw,
        'recordsTotal': records_total,
        'recordsFiltered': records_filtered,
        'data': data
    })
