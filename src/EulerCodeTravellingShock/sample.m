N=[10, 20,30,40,50,60];

error=[0.0386, 0.0315, 0.0247, 0.0202, 0.0170, 0.0149];

plot(log(N),log(error))

(log(error(6))-log(error(5)))/(log(1/60)-log(1/50));