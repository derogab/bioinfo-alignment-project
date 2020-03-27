# bioinfo-alignment-project
### Program to print reads from a SAM/BAM file filtered for start positions.

```
$ python app.py -h
usage: app.py [-h] -f FILE [-p POS [POS ...]] [--only-pos]

Print reads from a SAM/BAM file filtered for start positions.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  path to input file
  -p POS [POS ...], --pos POS [POS ...]
                        list of positions
  --only-pos            get list of all align start positions
```

##### W/ Docker 
```
docker build -t bioinfo-alignment-project .
```
```
docker run \
  -v /host/path/to/examples/:/usr/src/app/examples/ \
  bioinfo-alignment-project \
  [-h] -f FILE [-p POS [POS ...]] [--only-pos]
```

### Dependencies
###### Debian / Ubuntu
```
sudo apt-get install libbz2-dev liblzma-dev python3-dev samtools
```
```
sudo -H python3 -m pip install pysam
```