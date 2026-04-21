# Process helm charts playground

## What is a process ?

In order to have a broad compatibility, we chose to follow standard, and in the case of computation the [OGC API - process](https://ogcapi.ogc.org/processes/) corresponded to our needs. Their definition of a process is the following "The OGC API - Processes standard supports the wrapping of computational tasks into executable processes that can be offered by a server through a Web API and be invoked by a client application. The standard specifies a processing interface to communicate over a RESTful protocol using JavaScript Object Notation (JSON) encodings. Typically, these processes execute well-defined algorithms that ingest vector and/or coverage data to produce new datasets."

## Contribution steps

To add any process to our platform, you will have to make a merge request to this repositories.

> ⚠️ **This is a development repository**: A chart older than a week will be automatically deleted. To pass your process into production, please see [Create a merge request](#create-a-merge-request)

This tutorial provides guidelines to create your own helm chart in order to deploy your process on the EDITO datalab. Do not hesitate to look at [other charts](https://gitlab.mercator-ocean.fr/pub/edito-infra/process-helm-charts) to get inspired. In our example we use Python but any language can be used to create a process. Please adapt the python references to the other languages if needed.

First thing first, you will need to have a container image hosted on a public repository. This image should run a container exposing environment variables to determine inputs and outputs location, if they are needed. We strongly recommand to follow the [12-factor app](https://12factor.net/) methodology.

In this tutorial, we will take the [Coral bleaching detection process](https://datalab.dive.edito.eu/launcher/process-playground/coral-bleaching-detection) as an example. At the time being, it only has one input parameter determining if the process is launched as a small demonstration (small resources needed to run).

And as you can see, this project satisfies the minimal requirement for hosting a process on the datalab, which is having a public container image available with environment variables if an input is needed.

## Containerizing your computation code

In order to containerize our computation code, we used Docker.
If needed, please read the documentation to know more about getting started with [Docker](https://docs.docker.com/get-started/). In case you are using python together with a micromamba environment, [micromamba quick start](https://micromamba-docker.readthedocs.io/en/latest/quick_start.html) can be of help.

As a further example, here is one of our Dockerfile using python with a micromamba environment:

```Dockerfile
FROM mambaorg/micromamba:1.4.1-kinetic

COPY --chown=$MAMBA_USER:$MAMBA_USER ./conda_environment_bleaching.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

COPY ./bleaching.py /bleaching.py

CMD [ "python", "/bleaching.py" ]
```

The important points to follow is that you will need to reproduce the steps you are doing manually inside the Dockerfile in order for your code to be able to run. So all dependencies must be installed, via an environment manager, a dependency manager or one by one.
Furthermore the input must be made available via an environment variable, either its whole content or an URL to access it, and used accordingly into your code.
At last you must launch your code, as done in the example with `CMD [ "python", "/bleaching.py" ]`.

### Output written to user storage

The process used as a [template](https://gitlab.mercator-ocean.fr/pub/edito-infra/process-playground/-/tree/main/coral-bleaching?ref_type=heads) directly writes the output data to the user personal storage.
The code inside the container needs to write its output in a specific path determined thanks to an environment variable `EDITO_INFRA_OUTPUT`.
Thanks to the use of this environment variable the template process contains an additional step that will copy the content of this path to the user personal storage. Please go to [Writing directly to the user personal storage](#writing-directly-to-the-user-personal-storage) to know more about it.

### Python use of environment variable

In order to use the environment variables inside the python code, only the built-in `os` library is needed.
In this project the input is a boolean named `SMALL_DEMO` and used inside the python code to limit the quantity of data used for the demonstration.

```python
import os

os.environ.get("SMALL_DEMO")
```

## Clone the repository

Once you conternized your process, you can then clone the repository to start your contribution.

```sh
git clone https://gitlab.mercator-ocean.fr/pub/edito-infra/process-playground.git
```

## Create your own chart folder

You can start by copying the content of the `coral-bleaching` folder inside your own folder.

```sh
cp coral-bleaching my-process
```

If you know what you are doing, you can also start from scratch, with an empty helm chart.
```sh
helm create my-process
```

## Update the chart configuration

If you copied the `coral-bleaching` folder, you will then need to adjust some files.

### Edit the `Chart.yaml` file

Change the following fields and **leave the others unchanged**:

- **name** (the name of your process. This name must only consist of lower case alphanumeric characters, start with an alphabetic character, and end with an alphanumeric character. Hyphens (-) are allowed, but are known to be a little trickier to work with in Helm templates. The directory that contains a chart MUST have the same name as the chart)
- **description** (a brief description of your process)
- **home** (a page to learn more about your process, generate a "Learn more" button on the process tile)
- **icon** (an image that represent the underlying process)
- **keywords** (a list of useful keywords that can be used to retrieve your process from the datalab search bar)
- **version** (the version of the chart. Starts with 1.0.0 and update later if you need some changes)
- **appVersion** (the version of the process running inside your docker container. Maybe a version of your computation is present inside the repository where your process is versioned)

**All of these attributes are mandatory**, please find an icon even a generic one to illustrate your process.

### Edit the `templates/NOTES.txt` file

The content will be rendered and displayed in a pop-up window while the process is being launched. This text targets any user discovering your process. If you have an estimation of the time needed for the process to complete, it could be an interesting information to add.
By default the name of your process will be indicated if you keep the original notes, you can have access to Helm values such as `{{ .Chart }}`, `{{ .Release }}`, etc., as you can see in our example.
You can use other Helm values in this template file. Please take a look at the official [Helm documentation](https://helm.sh/docs/chart_template_guide/notes_files/) to learn more about it.

### Edit the `values.yaml` file

Replace the input environment variable `smallDemo` by your own and its default value, or add new ones. If you need an output environment variable, add a "processOutputs" section at the bottom with its name. As you can see, this named variable is not exactly the same one as the name from the python code. The correspondence will be done in a further step.

```yaml
...
demo:
  smallDemo: true
...
```

### Edit the `values.schema.json` file

Replace into the file the input information by your own and add any necessary according to the environment variables you put into `value.yaml`.

```json
{
    ...
        "demo": {
          "description": "Process inputs",
          "type": "object",
          "properties": {
            "smallDemo": {
              "type": "boolean",
              "description": "To run a small demo of the process",
              "default": true
            }
          }
        }
    ...
}
```

### Edit the `templates/job.yaml` file

Replace the container image reference `image: docker.mercator-ocean.fr/moi-docker/bleaching:...` by your own. Replace or add the environment variables by precising their names and values. The easiest way is to use the Helm values to directly inject the value of the environment variable in this way : `"{{ .Values.demo.smallDemo }}"`.
As you can see, the link between the environment variables and the names used inside the `values.yaml` file is indicated here.

 ```yaml
        env:
        - name: SMALL_DEMO
          value: "{{ .Values.demo.smallDemo }}"
 ```
When you push your branch, your charts will automatically be published and accessible on EDITO datalab [process playground](https://datalab.dive.edito.eu/process-catalog/process-playground) (there may be a 5-minute refresh delay).

**Note**: in the `job.yaml` file, the following lines must not be touched:
```yaml
metadata:
  name: {{ .Release.Name }}
```

## Create a merge request

Once you think your chart is ready to be published, you can:

1. Make sure the metadata are complete in the `Chart.yaml` and `README.md` files.
2. Please provide somehow a point of contact for the users to reach you.
3. Pick a catalog category in which your contribution fit the best.
4. Create a merge request on the repository and ping @pub/edito-infra/codeowners in the description to catch our attention.

If everything is good, we will migrate your chart to the category you provided, and you will be granted accesses to maintain them (bug fixes, new versions, etc.).

## Additional information

### Writing directly to the user personal storage

Thanks to the use of `EDITO_INFRA_OUTPUT` environment variable the template process contains an additional step that will copy the content of this path to the user personal storage. This step can be found in the file [job.yaml](https://gitlab.mercator-ocean.fr/pub/edito-infra/process-playground/-/blob/main/coral-bleaching/templates/job.yaml?ref_type=heads) as a container named `copy-output`.

These variables are added as environment variables in the container with:
```yaml
envFrom:
{{- if .Values.s3.enabled }}
- secretRef:
    name: {{ include "library-chart.secretNameS3" . }}
{{- end }}
```

This container needs to have access to specific environment variables to be able to access the user S3 bucket:

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `AWS_SESSION_TOKEN`
- `AWS_S3_ENDPOINT`
- `AWS_DEFAULT_REGION`

These variables have been exported thanks to the `s3` section of the [values.schema.json](https://gitlab.mercator-ocean.fr/pub/edito-infra/process-playground/-/blob/main/coral-bleaching/values.schema.json?ref_type=heads) and the presence of the file [secret-s3.yaml](https://gitlab.mercator-ocean.fr/pub/edito-infra/process-playground/-/blob/main/coral-bleaching/templates/secret-s3.yaml?ref_type=heads).


### Include Copernicus Marine Service credentials
It is possible to load Copernicus Marine Service credentials as environement variables in the process. The following configuration will automatically import the credentials configured in the user's `My Account`.

First, to automatically load Copernicus Marine Service credentials into the service configuration, add the following property in the `values.schema.json` file:
```json
{
  "properties": {
    ...
    "copernicusMarine": {
      "x-onyxia": {
        "overwriteSchemaWith": "copernicusMarine.json"
      }
    },
    ...
  }
}
```

Add the following properties in the `values.yaml` files:
```yaml
copernicusMarine:
  enabled: false
  username: ""
  password: ""
```

Then create a `secret-copernicusmarine.yaml` file inside the `templates` folder with the following content:
```yaml
{{- define "library-chart.secretNameCopernicusMarine" -}}
{{- if .Values.copernicusMarine.enabled }}
{{- $name:= (printf "%s-secretcopernicusmarine" (include "library-chart.fullname" .) )  }}
{{- default $name .Values.copernicusMarine.secretName }}
{{- else }}
{{- default "default" .Values.copernicusMarine.secretName }}
{{- end }}
{{- end }}

{{- if .Values.copernicusMarine.enabled -}}
apiVersion: v1
kind: Secret
metadata:
  name: {{ include "library-chart.secretNameCopernicusMarine" . }}
  labels:
    {{- include "library-chart.labels" . | nindent 4 }}
stringData:
  COPERNICUSMARINE_SERVICE_USERNAME: "{{ .Values.copernicusMarine.username }}"
  COPERNICUSMARINE_SERVICE_PASSWORD: "{{ .Values.copernicusMarine.password }}"
{{- end }}
```

Finally, load the secret values as environment variables in the container:
```yaml
envFrom:
  {{- if .Values.copernicusMarine.enabled }}
  - secretRef:
      name: {{ include "library-chart.secretNameCopernicusMarine" . }}
  {{- end }}
```

### GPU-based processes

If the process relies on a GPU card, the container must include a cuda integration, for example it could derive a micromamba image integrating cuda:
```
# Dockerfile
FROM mambaorg/micromamba:1.5.6-focal-cuda-12.1.1
```

In addition, some changes must be applied to `values.schema.json` and `template/job.yaml`:

In `values.schema.json`, use `ide/resources-gpu.json` instead of `ide/resources.json`:
```
...
  "properties": {
    "resources": {
      "x-onyxia": {
        "overwriteSchemaWith": "ide/resources-gpu.json"
      }
    },
    ...
```

In `template/job.yaml`, add the following `tolerations`:
```
...
spec:
  template:
    spec:
      tolerations:
      - effect: NoSchedule
        key: node.cloudferro.com/type
        operator: Equal
        value: gpu
    ...
```
