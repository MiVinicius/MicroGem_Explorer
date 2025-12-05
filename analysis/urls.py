from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('analyses/', views.analysis_list, name='analysis_list'),
    path('create/', views.analysis_create, name='analysis_create'),
    path('<int:pk>/', views.analysis_detail, name='analysis_detail'),
    path('<int:pk>/delete/', views.analysis_delete, name='analysis_delete'),
    path('<int:pk>/variants_api/', views.analysis_variants_api, name='analysis_variants_api'),
]
