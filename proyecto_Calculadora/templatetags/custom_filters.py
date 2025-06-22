from django import template
import re

register = template.Library()

@register.filter(name='replace_bold')
def replace_bold(value):
    """Reemplaza **texto** con <strong>texto</strong>."""
    return re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', value)