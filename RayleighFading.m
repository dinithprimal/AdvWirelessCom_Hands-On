function h = RayleighFading(Nr);

   x = sqrt(1/2).*rand(1,Nr);
   y = sqrt(1/2).*rand(1,Nr);
   
   h = abs(x + i*y);
   
end