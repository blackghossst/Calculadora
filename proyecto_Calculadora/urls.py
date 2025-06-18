from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('metodo_muller', views.resolver, name='resolver'),
    path('formulario_muller', views.mostrar_formulario_muller, name='formulario_muller'),
]
