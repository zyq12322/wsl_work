g++ -std=c++17 -O3 -march=native     -Wall -Wextra     var_op_function.cpp     dealiaser.cpp     DNO_calc_OE.cpp     DNO_output_control.cpp     DNO_time_integral.cpp     DNO_time_integral_main.cpp     -o dno_sim     -lfftw3 -lm
./dno_sim
python3 dno_animate.py --out result.mp4 --eta-scale 50 --beta-offset 10