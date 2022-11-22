from django.shortcuts import render
from django.http import HttpResponse
from .forms import UploadFileForm

# Create your views here.
# def index(response):
#     return render(response, "main/base.html", {})

def upload_file(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        for f in request.FILES.getlist('file'):
            print(str(f)) 
        file = request.FILES.getlist('file')
        return HttpResponse("The name of this file is " + str(file))
    else:
        form = UploadFileForm()
    return render(request, "main/base.html", {'form': form})

def v1(response):
    return HttpResponse("<h1>View 1</h1>")