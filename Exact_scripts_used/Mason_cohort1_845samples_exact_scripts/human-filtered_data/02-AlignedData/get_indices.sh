echo "get indices"

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME


fqPath=/path/to/01-TrimmedData/Fastq_RG


for sample in $(cat samples_H5MK2DSXY_1-2.txt)
do
        cat $fqPath/${sample}/${sample}_H5MK2DSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5MK2DSXY_1-2.txt
done


for sample in $(cat samples_H5MK2DSXY_3-4.txt)
do
        cat $fqPath/${sample}/${sample}_H5MK2DSXY_3_R1.report.json | grep '"multiplex_barcode":' >> indices_H5MK2DSXY_3-4.txt
done


for sample in $(cat samples_H5MNLDSXY_1-2.txt)
do
        cat $fqPath/${sample}/${sample}_H5MNLDSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5MNLDSXY_1-2.txt
done


for sample in $(cat samples_H5MNLDSXY_3-4.txt)
do
        cat $fqPath/${sample}/${sample}_H5MNLDSXY_3_R1.report.json | grep '"multiplex_barcode":' >> indices_H5MNLDSXY_3-4.txt
done


for sample in $(cat samples_H5NJ5DSXY_1-2.txt)
do
        cat $fqPath/${sample}/${sample}_H5NJ5DSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NJ5DSXY_1-2.txt
done


for sample in $(cat samples_H5NJ5DSXY_3.txt)
do
        cat $fqPath/${sample}/${sample}_H5NJ5DSXY_3_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NJ5DSXY_3.txt
done


for sample in $(cat samples_H5NKGDSXY_1-2.txt)
do
	cat $fqPath/${sample}/${sample}_H5NKGDSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NKGDSXY_1-2.txt
done


for sample in $(cat samples_H5NKGDSXY_3-4.txt)
do
        cat $fqPath/${sample}/${sample}_H5NKGDSXY_3_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NKGDSXY_3-4.txt
done


for sample in $(cat samples_H5NTTDSXY_1-2.txt)
do
        cat $fqPath/${sample}/${sample}_H5NTTDSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NTTDSXY_1-2.txt
done


for sample in $(cat samples_H5NTTDSXY_3-4.txt)
do
        cat $fqPath/${sample}/${sample}_H5NTTDSXY_3_R1.report.json | grep '"multiplex_barcode":' >> indices_H5NTTDSXY_3-4.txt
done


for sample in $(cat samples_H5YJHDSXY_1-2.txt)
do
        cat $fqPath/${sample}/${sample}_H5YJHDSXY_1_R1.report.json | grep '"multiplex_barcode":' >> indices_H5YJHDSXY_1-2.txt
done


end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"
