
# How to create a new template for cookiecutter

+ 1. Get started, install `Cookiecutter` using `pip` or `conda`

```
pip install cookiecutter
```

+ 2. Create a new directory for your template


```
mkdir new_project
cd new_project
```

+ 3. Create template structures

```
mkdir {{cookiecutter.project_name}}
cd {{cookiecutter.project_name}}
mkdir {{cookiecutter.documentation}}
mkdir {{cookiecutter.data}}
mkdir {{cookiecutter.source}}
mkdir {{cookiecutter.bin}}
mkdir {{cookiecutter.results}}
cd {{cookiecutter.data}}
mkdir {{cookiecutter.raw_data}}
mkdir {{cookiecutter.clean_data}}
mkdir {{cookiecutter.sample_info}}
cd ..
```

+ 4. Create a README file inside each project directory

```
wget https://gist.githubusercontent.com/PurpleBooth/109311bb0361f32d87a2/raw/824da51d0763e6855c338cc8107b2ff890e7dd43/README-Template.md -O tmp.md
cat tmp.md | sed 's/Project Title/{{cookiecutter.project_name}}/' > {{cookiecutter.README}}.md
rm -f tmp.md
```

+ 5. Finally, we need to create the `cookiecutter.json` file, like this:

```
{
    "project_name": "new_project",
    "documentation": "doc",
    "data": "data",
    "clean_dta": "clean_data",
    "raw_data": "raw_data",
    "sample_info": "sample_info",
    "source": "src",
    "bin": "bin",
    "results": "results",
    "README": "README"
}
```

If you followed the steps above, you should have this directory structure inside the `new_project` directory.

```
.
|-- cookiecutter.json
|-- {{cookiecutter.project_name}}
|   |-- {{cookiecutter.bin}}
|   |-- {{cookiecutter.data}}
|   |   |-- {{cookiecutter.clean_data}}
|   |   |-- {{cookiecutter.raw_data}}
|   |   `-- {{cookiecutter.sample_info}}
|   |-- {{cookiecutter.documentation}}
|   |-- {{cookiecutter.README}}.md
|   |-- {{cookiecutter.results}}
|   `-- {{cookiecutter.source}}
```
+ 6. Now we can create new projects using this template, suppose the path of `new_project` is `~/biodata/new_project/`, siwtch to the directory where you want to create the a project.

```
cookiecutter ~/biodata/new_project
project_name [new_project]: msms
documentation [doc]: 
data [data]: 
raw_data [raw_data]: 
clean_data [clean_data]: 
sample_info [sample_info]: 
source [src]: 
bin [bin]: 
results [results]: 
README [README]: 
```
+ 7. Now your folder should looks like this:

```
temp001/
|-- bin
|-- data
|   |-- clean_data
|   |-- raw_data
|   `-- sampel_info
|-- doc
|-- README.md
|-- results
`-- src

8 directories, 1 file
```



