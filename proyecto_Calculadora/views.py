from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import render
import numpy as np
import json
from django.db import IntegrityError
import re
from sympy import sympify, Symbol, lambdify, E
from .forms import AdaptiveGaussianQuadratureForm


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

from django.shortcuts import render
from .forms import AdaptiveGaussianQuadratureForm
from sympy import sympify, Symbol, lambdify, E, pi
import numpy as np

def adaptive_gaussian_quadrature_calculator(request):
    result = None
    error = None
    steps = []

    if request.method == 'POST':
        form = AdaptiveGaussianQuadratureForm(request.POST)
        if form.is_valid():
            func_str = form.cleaned_data['function']
            a = form.cleaned_data['a']
            b = form.cleaned_data['b']

            try:
                if a >= b:
                    raise ValueError("El límite inferior debe ser menor que el límite superior.")

                x = Symbol('x')
                expr = sympify(func_str, locals={"e": E, "pi": pi})
                f = lambdify(x, expr, modules=['numpy'])

                steps.append(f"**Función:** $f(x) = {func_str}$")
                steps.append(f"**Intervalo de integración:** $[{a}, {b}]$")
                steps.append("---") # Separador visual

                def gauss_quad_2pts(f, a, b, profundidad):
                    # No generamos 'pasos' internos para simplificar la salida
                    puntos = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
                    pesos = np.array([1, 1])

                    mitad = (a + b) / 2
                    longitud = (b - a) / 2

                    suma = 0
                    x_evals_str = []
                    fx_evals_str = []
                    for i, (p, w) in enumerate(zip(puntos, pesos)):
                        x_eval = mitad + longitud * p
                        fx = f(x_eval)
                        x_evals_str.append(f"x_{i+1} = {x_eval:.6f}")
                        fx_evals_str.append(f"f(x_{i+1}) = f({x_eval:.6f}) = {fx:.6f}")
                        suma += w * fx

                    resultado = longitud * suma
                    # Retornamos los detalles para que la función llamadora los imprima si es necesario
                    return resultado, x_evals_str, fx_evals_str


                def metodo_adaptativo(f, a, b, pasos, profundidad=0, max_prof=10):
                    indent = '  ' * profundidad # Indentación para la visualización de la profundidad
                    
                    pasos.append(f"{indent}**Analizando subintervalo:** $[{a:.4f}, {b:.4f}]$ (Profundidad: {profundidad})")

                    # Cálculo de I1 (integral del intervalo completo con 2 puntos)
                    I1, x_I1, fx_I1 = gauss_quad_2pts(f, a, b, profundidad)
                    pasos.append(f"{indent}  * Cálculo inicial con 2 puntos (I₁):")
                    pasos.append(f"{indent}     Puntos: {', '.join(x_I1)}")
                    pasos.append(f"{indent}     Valores de la función: {', '.join(fx_I1)}")
                    pasos.append(f"{indent}     Integral I₁ ≈ {I1:.6f}")

                    medio = (a + b) / 2
                    
                    # Cálculo de I2 (suma de integrales de los dos subintervalos)
                    I2_izq, x_I2_izq, fx_I2_izq = gauss_quad_2pts(f, a, medio, profundidad + 1)
                    I2_der, x_I2_der, fx_I2_der = gauss_quad_2pts(f, medio, b, profundidad + 1)
                    I2 = I2_izq + I2_der

                    pasos.append(f"{indent}  * Dividimos en $[{a:.4f}, {medio:.4f}]$ e $[{medio:.4f}, {b:.4f}]$:")
                    pasos.append(f"{indent}     Integral izquierda (I_izq) ≈ {I2_izq:.6f}")
                    pasos.append(f"{indent}     Integral derecha (I_der) ≈ {I2_der:.6f}")
                    pasos.append(f"{indent}     Integral I₂ = I_izq + I_der ≈ {I2_izq:.6f} + {I2_der:.6f} = {I2:.6f}")

                    error_aprox = abs(I2 - I1)
                    pasos.append(f"{indent}  **Error estimado:** $|I_2 - I_1| = |{I2:.6f} - {I1:.6f}| = {error_aprox:.6f}")

                    # Tolerancia y profundidad máxima (pueden ser parámetros de entrada del formulario)
                    TOLERANCIA = 1e-4 
                    MAX_PROFUNDIDAD = 10

                    if error_aprox < TOLERANCIA or profundidad >= MAX_PROFUNDIDAD:
                        pasos.append(f"{indent}  ✔️ **Aceptado:** El error es menor que la tolerancia ({TOLERANCIA}) o se alcanzó la profundidad máxima.")
                        pasos.append(f"{indent}     Resultado para este subintervalo: ≈ {I2:.6f}")
                        pasos.append(f"{indent}---") # Separador
                        return I2
                    else:
                        pasos.append(f"{indent}  ⚠️ **Subdividiendo:** El error ({error_aprox:.6f}) es mayor que la tolerancia ({TOLERANCIA}).")
                        pasos.append(f"{indent}---") # Separador
                        izq = metodo_adaptativo(f, a, medio, pasos, profundidad + 1)
                        der = metodo_adaptativo(f, medio, b, pasos, profundidad + 1)
                        return izq + der

                # Ejecutar el cálculo
                integral = metodo_adaptativo(f, a, b, steps)
                result = round(integral, 6) # Mayor precisión en el resultado final
                steps.append(f"---")
                steps.append(f"## ✅ **Resultado Final de la Integral:** $\int_{a}^{b} {func_str} \,dx \\approx {result}$")

            except Exception as e:
                error = f"❌ **Error durante el cálculo:** {str(e)}"
        else:
            error = "Por favor, corrige los errores en el formulario."
    else:
        form = AdaptiveGaussianQuadratureForm()

    return render(request, 'cuadratica.html', {
        'form': form,
        'result': result,
        'error': error,
        'steps': steps
    })