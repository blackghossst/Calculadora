from django.db import models


class proyecto_Calculadora(models.Model):
    nombre = models.CharField(max_length=100)
    tipo = models.CharField(max_length=50)

class IntegralRegistro(models.Model):
    funcion = models.CharField(max_length=200)
    limite_inferior = models.FloatField()
    limite_superior = models.FloatField()
    tolerancia = models.FloatField(default=1e-6)
    resultado = models.FloatField()
    fecha = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"∫({self.funcion}) dx [{self.limite_inferior}, {self.limite_superior}] ≈ {self.resultado}"
