from django import forms


class IGVform(forms.Form):
    sample_name = forms.CharField(max_length=100)
    run_name = forms.CharField(max_length=100)
    reference = forms.CharField(max_length=100)
    unique_id = forms.CharField(max_length=100)
    project_name = forms.CharField(max_length=100)


class download_form(forms.Form):
    file_path = forms.CharField(max_length=300)

    class Meta:

        widgets = {
            "myfield": forms.TextInput(
                attrs={"style": "border-color:darkgoldenrod; border-radius: 10px;"}
            ),
        }
