from django.urls import path

from . import views
from .views import upload_file

urlpatterns = [
    path('', views.v1, name='index'),
    path('v1/', views.v1, name='view 1'),
    path('list/', upload_file, name='list')
]