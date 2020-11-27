set term pdfcairo monochrome font "dejaVu Sans Mono,10"
set out "OUTPUT.pdf"
set xlabel "positions (x) in 1D Box"
set ylabel "wave-functions"
set title "Wave Functions for Particle in a Double Well Problem"
set xtics nomirror  
set ytics nomirror  
set key top right
#set colorsequence [default | classic |podo] 
#f(x)=a*(1-exp(-b*x))
#f(x)= a*exp(b*x)- c*exp(-d*x)+e
#b = 0.5
#a = 50
#c = 50
#d = 0.5
#e = 20
#fit f(x) "/home/dipankar/heat.txt" u 1:2 via a,b,c,d,e
#f1(x)= 73.4523*exp(0.1*x)- 53.4523*exp(-0.1*x)+20
#ff(x)=4*(exp(0.8*x)-exp(-0.5*x))/1.3+2*exp(-0.5*x)
#plot f(x) w l lt 2 t sprintf("f(x)= %.4f*exp(%.4f*x)- %.4f*exp(-%.4f*x)+%.4f",a,b,c,d,e), "/home/dipankar/heat.txt" u 1:2 w p ls 1 notitle
V(x)=b0+b1*(x-x0)+b2*(x-x0)*(x-x1)+b3*(x-x0)*(x-x1)*(x-x2)+b4*(x-x0)*(x-x1)*(x-x2)*(x-x3)
 b0=0
 b1=-5
 b2=5
 b3=-3
 b4=1.3333
 x0=-2
 x1=-1
 x2=0
 x3=1


plot "/home/dipankar/heat.txt" u 1:2 w l lt 1 notitle, "/home/dipankar/heat.txt" u 1:2 w p pt 5 ps 0.2 notitle, V(x) w l lt 2 t "V(x)"
