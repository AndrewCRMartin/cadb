for file in *.c *.h Makefile 00READ.ME
do
   ci -l $file
done
