from django import forms

class AdaptiveGaussianQuadratureForm(forms.Form):
    function = forms.CharField(
        label="Función f(x)",
        initial="x**2",
        help_text="Introduce la función en términos de 'x' (ej: x**2 + sin(x), exp(x), 1/x)."
    )
    a = forms.FloatField(
        label="Límite inferior (a)",
        initial=0.0,
        help_text="Valor inicial del intervalo de integración."
    )
    b = forms.FloatField(
        label="Límite superior (b)",
        initial=1.0,
        help_text="Valor final del intervalo de integración."
    )