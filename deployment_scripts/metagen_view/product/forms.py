from django import forms


class UploadFileForm(forms.Form):
    technology = forms.ChoiceField(
        choices=[("illumina", "Illumina"), ("nanopore", "nanopore")]
    )
    file_r1 = forms.FileField(
        widget=forms.ClearableFileInput(attrs={"multiple": True, "accept": "*"})
    )
    file_r2 = forms.FileField(
        required=False, widget=forms.FileInput(attrs={"accept": "application/fastq"})
    )
    name = forms.CharField(
        max_length=254, widget=forms.TextInput(attrs={"placeholder": "Enter your name"})
    )

    project_name = forms.CharField(
        max_length=254,
        widget=forms.TextInput(attrs={"placeholder": "Enter your project name"}),
    )
