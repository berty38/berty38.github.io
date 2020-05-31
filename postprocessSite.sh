cat index.html | sed 's/<\/head>.*//' > tmp1.txt

cat index.html | sed 's/.*<\/head>/<\/head>/' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > index.html


cat Bert_Huang/Home.html | sed -n '1,/<\/head>/p' > tmp1.txt

cat Bert_Huang/Home.html | sed -n '1,/<\/head>/!p' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > Bert_Huang/Home.html


cat Bert_Huang/Research.html | sed -n '1,/<\/head>/p' > tmp1.txt

cat Bert_Huang/Research.html | sed -n '1,/<\/head>/!p' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > Bert_Huang/Research.html


cat Bert_Huang/Teaching.html | sed -n '1,/<\/head>/p' > tmp1.txt

cat Bert_Huang/Teaching.html | sed -n '1,/<\/head>/!p' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > Bert_Huang/Teaching.html


cat Bert_Huang/Hobbies.html | sed -n '1,/<\/head>/p' > tmp1.txt

cat Bert_Huang/Hobbies.html | sed -n '1,/<\/head>/!p' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > Bert_Huang/Hobbies.html


cat Bert_Huang/Links.html | sed -n '1,/<\/head>/p' > tmp1.txt

cat Bert_Huang/Links.html | sed -n '1,/<\/head>/!p' > tmp2.txt

cat tmp1.txt analytics.txt tmp2.txt > Bert_Huang/Links.html

python addWordpress.py > tmp1.txt
mv tmp1.txt Bert_Huang/Blog.html
