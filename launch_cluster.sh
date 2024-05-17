# for i in {0..4}
# do
#     koef=$(echo "0.75+$i*0.05" | bc)
#     python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0$koef.json --outputFolder pion_100k_300GeV_mc_Ekoef-0$koef-abs1-nop
# done

# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.75.json --outputFolder pion_100k_300GeV_mc_Ekoef-075
# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.75.json --outputFolder pion_100k_300GeV_mc_Ekoef-075
# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.75.json --outputFolder pion_100k_300GeV_mc_Ekoef-075
# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.75.json --outputFolder pion_100k_300GeV_mc_Ekoef-075
# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.75.json --outputFolder pion_100k_300GeV_mc_Ekoef-075

# python run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.50.json --outputFolder pion_100k_300GeV_mc_Ekoef-0.50

python3 run_submissions_template.py --nEvents 100000 --jobs 100 --param arguments_0.90_180.json --outputFolder pion_100k_180GeV_mc_Ekoef-0.90-abs1-nop