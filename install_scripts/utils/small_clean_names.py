import sys

fq1 = sys.argv[1]
f2 = sys.argv[2]


def clean_fastq_headers_python(fastq, temp_fq):
    """
    Clean fastq header using python.
    """
    counter = 0  # counter to keep track of headers
    with open(temp_fq, "w") as f:
        for line in open(fastq):
            line = line.strip()
            if line.startswith("@") and counter == 0:
                if line[-2:] == "/1" or line[-2:] == "/2":
                    line = line[:-2]

            f.write(line + "\n")

            counter += 1
            if counter == 4:
                counter = 0


clean_fastq_headers_python(fq1, f2)
