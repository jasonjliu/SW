# Smith Waterman Local Alignment

This algorithm implementation performs the Smith Waterman Local Alignment with Open and Extension Gaps.

The `blosum62.txt` file provides the scores for match and mismatch.

The usage is as follows:

```
python hw1.py -i <input file> -s <score file> > output.txt
```

Additionally, open Gap (`openGap`) and extension Gap (`extGap`) scores can be defined as follows:

```
python hw1.py -i <input file> -s <score file> -o <openGap> -e <extGap> > output.txt

```

The default values for `openGap` and `extGap` are `-2` and `-1`, respectively.

The first example is with the input file `sample-input1.txt`

It should produce `sample-output1.txt`
```
# run the algorithm for sample-input1.txt and produce the file output1.txt
python hw1.py -i sample-input1.txt -s blosum62.txt > output1.txt
# the files could be compared as
diff -E -b output1.txt sample-output1.txt
```

The second example is with the input file `sample-input2.txt`

It should produce `sample-output2.txt`
```
# run the algorithm for sample-input2.txt and produce the file output2.txt
python hw1.py -i sample-input2.txt -s blosum62.txt > output2.txt
# the files could be compared as
diff -E -b output2.txt sample-output2.txt
```

Finally, you could test the algorithm on a new input file, `input.txt`, as follows:
```
python hw1.py -i input.txt -s blosum62.txt > output_new.txt

```

