from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import render
import numpy as np
import json
from django.db import IntegrityError
import re

def home(request):
    return render(request, 'home.html')

def corregir_ecuacion(equation):
    equation = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', equation)
    equation = re.sub(r'([a-zA-Z])\(', r'\1*(', equation)
    equation = re.sub(r'\)([a-zA-Z])', r')*\1', equation)
    equation = re.sub(r'\)\(', r')*(', equation)
    return equation

def formato_respuesta(numero):
    return f"{numero:.4f}".rstrip('0').rstrip('.') if '.' in f"{numero:.4f}" else f"{numero}"


def resolver(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        ecuacion = data.get('equation')
        ecuacion = corregir_ecuacion(ecuacion)
        x0 = float(data.get('x0'))
        x1 = float(data.get('x1'))
        x2 = float(data.get('x2'))

        try:
            def funcion(x):
                return eval(ecuacion, {"x": x, "np": np, "__builtins__": {}})

            resultado = metodo_muller(funcion, x0, x1, x2)

            if 'error' in resultado:
                return JsonResponse({'error': resultado['error']})
            else:
                return JsonResponse({
                    'raiz': resultado['raiz'],
                    'iteraciones': resultado['iteraciones']
                })

        except Exception as e:
            return JsonResponse({'error':'Error al procesar la ecuación: ingrese solo numeros y operadores matemáticos.'})

    return JsonResponse({'error': 'Método no permitido'})

def metodo_muller(f, x0, x1, x2, tol=1e-6, max_iter=100):
    iteraciones = []

    for i in range(max_iter):
        h0 = x1 - x0
        h1 = x2 - x1
        sigma0 = (f(x1) - f(x0)) / h0
        sigma1 = (f(x2) - f(x1)) / h1
        a = (sigma1 - sigma0) / (h1 + h0)
        b = a * h1 + sigma1
        c = f(x2)
        discriminante = b ** 2 - 4 * a * c

        if discriminante < 0:
            return {"error": "Discriminante negativo. No se puede encontrar una raíz real."}

        D = np.sqrt(discriminante)

        if abs(b + D) > abs(b - D):
            E = b + D
        else:
            E = b - D

        if E == 0:
            return {"error": "División por cero. Intenta con otros puntos iniciales."}

        h = -2 * c / E
        x3 = x2 + h

        iteracion_info = {
            'iteracion': i + 1,
            'h0': f"x₁ - x₀ = {formato_respuesta(x1)} - {formato_respuesta(x0)} = {formato_respuesta(h0)}",
            'h1': f"x₂ - x₁ = {formato_respuesta(x2)} - {formato_respuesta(x1)} = {formato_respuesta(h1)}",
            'sigma0': f"(f(x₁) - f(x₀)) / h₀ = ({formato_respuesta(f(x1))} - {formato_respuesta(f(x0))}) / {formato_respuesta(h0)} = {formato_respuesta(sigma0)}",
            'sigma1': f"(f(x₂) - f(x₁)) / h₁ = ({formato_respuesta(f(x2))} - {formato_respuesta(f(x1))}) / {formato_respuesta(h1)} = {formato_respuesta(sigma1)}",
            'a': f"(σ₁ - σ₀) / (h₁ + h₀) = ({formato_respuesta(sigma1)} - {formato_respuesta(sigma0)}) / ({formato_respuesta(h1)} + {formato_respuesta(h0)}) = {formato_respuesta(a)}",
            'b': f"a * h₁ + σ₁ = {formato_respuesta(a)} * {formato_respuesta(h1)} + {formato_respuesta(sigma1)} = {formato_respuesta(b)}",
            'c': f"f(x₂) = {formato_respuesta(c)}",
            'discriminante': f"b² - 4ac = {formato_respuesta(b)}² - 4 * {formato_respuesta(a)} * {formato_respuesta(c)} = {formato_respuesta(discriminante)}",
            'D': f"√Discriminante = √{formato_respuesta(discriminante)} = {formato_respuesta(D)}",
            'E': f"b + D = {formato_respuesta(E)}",
            'h': f" -2 * c / E = -2 * {formato_respuesta(c)} / {formato_respuesta(E)} = {formato_respuesta(h)}",
            'xi_1': f" x₂ + h = {formato_respuesta(x2)} + {formato_respuesta(h)} = {formato_respuesta(x3)}",
            'f(xi_1)': f" {formato_respuesta(f(x3))}"
        }
        iteraciones.append(iteracion_info)

        if abs(h) < tol:
            return {"raiz": formato_respuesta(x3), "iteraciones": iteraciones}

        x0, x1, x2 = x1, x2, x3

    return {"error": "No se encontró la raíz en el número de iteraciones especificado.", "iteraciones": iteraciones}



def mostrar_formulario_muller(request):
    return render(request, 'metodo_muller.html')