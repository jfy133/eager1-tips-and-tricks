# How to debug in EAGER?
## Debugging setup
### EAGER-GUI won't generate the output folders
You have finishing setting up all the parameters, you press *Generate Config File* and either EAGER-GUI throws you an error or no error is thrown but when you check the output directory there is no output folder generated. Usually this EAGER behavious is due to an erroneous setup.

You can check a few things:

**Make sure that you have chosen an the right output folder:** EAGER will generate the output in whichever output folder you have indicated. If you have selected the another folder as output folder than the one you want the ouput to be, you will find there the output folder for EAGER. In this case, you can either:
1. Do the set up again selecting the correct output folder
2. Move the folders to the correct folder and change the EAGER config file (file ending with *.xml*). 
  1. You can do this manually in a text editor, changing these tags to: 
```
**<GUI__filepathresults>/path/to/wrong/output_folder/</GUI__filepathresults>** to **<GUI__filepathresults>/path/to/correct/output_folder/</GUI__filepathresults>**

**<GUI__resultspath>/path/to/wrong/output_folder/Sample</GUI__resultspath>** to **<GUI__resultspath>/path/to/correct/output_folder/Sample</GUI__resultspath>**
```
  2. Or you can replace the wrong path using *sed* in the following manner: 
     ```
     sed -i 's$<GUI__filepathresults>/path/to/wrong/output_folder/$<GUI__filepathresults>/path/to/correct/output_folder/$'
     sed -i 's$<GUI__resultspath>/path/to/wrong/output_folder/Sample$<GUI__resultspath>/path/to/correct/output_folder/Sample/$'
     ```

**Make sure you selected the correct type of data**
You can also experience the following: you set up EAGER-GUI, your press *Generate Config File*, no error is thrown, but when you check your output folder only one folder has been created. This can happen if you have **all** the *.fastq* files in a single folder. EAGER-GUI expects to have the input stored in one folder per sample and in each of those the corresponding *.fastq* like this

```
Input folder/
  Sample1/
    Sample1_L001_R1_001.fastq.gz
    Sample1_L001_R2_001.fastq.gz
  Sample2/
    Sample2_L001_R1_001.fastq.gz
    Sample2_L001_R2_001.fastq.gz
    .
    .
    .
  SampleN/
    SampleN_L001_R1_001.fastq.gz
    SampleN_L001_R2_001.fastq.gz
```

### EAGER-GUI only generates one folder

## Debugging EAGER-CLI
### EAGER runs stop early
#### X11 turn on
#### VCF2Genome is ticked on

## Results are not what I expected
### Check configuration parameters
