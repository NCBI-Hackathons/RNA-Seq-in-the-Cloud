gcloud container clusters create study-gtf --zone=us-east1-b --enable-autoscaling --preemptible --min-nodes=1 --max-nodes=100 --num-nodes=1 --scopes storage-rw --disk-type=pd-ssd --disk-size=200 --machine-type=n1-highcpu-16
gcloud container clusters get-credentials study-gtf --zone=us-east1-b

tail -n +2 stringtie_merge_samples.tsv | tr ' ' '_' | while IFS=$'\t' read -r study tissue runs; do
  name=$(echo -n "${study}-${tissue}-gtf" | tr [A-Z] [a-z]| tr -c [a-z0-9-] '-')
  echo $name

  if [ -f ${name}.yaml ] ; then continue; fi
  cat job_template.yaml | sed -r "s/\{name\}/$name/; s/\{study\}/$study/; s/\{tissue\}/$tissue/; s/\{runs\}/$runs/" > ${name}.yaml
  kubectl create -f ${name}.yaml
done
