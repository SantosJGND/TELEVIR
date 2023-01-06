from django import forms


class IGVform(forms.Form):
    project_name = forms.CharField(max_length=100)
    sample_name = forms.CharField(max_length=100)
    run_name = forms.CharField(max_length=100)
    reference = forms.CharField(max_length=100)
    unique_id = forms.CharField(max_length=100)


class download_form(forms.Form):
    file_path = forms.CharField(max_length=100)
