from django.db import models


class proyecto_Calculadora(models.Model):
    nombre = models.CharField(max_length=100)
    tipo = models.CharField(max_length=50)
