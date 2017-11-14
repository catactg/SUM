% ne fonctionne que pour des dimensions paires
function [fft3k_d fft3k] = noyau3D_PAIR(grille,lambda)
  longX = grille.long(1);
  longY = grille.long(2);
  longZ = grille.long(3);
  pas = grille.pas;
  k = zeros(grille.long);
  demilongX = longX/2;
  demilongY = longY/2;
  demilongZ = longZ/2;
  lambda2 = lambda*lambda;
  for i = 0:demilongX
   for j = 0:demilongY
    for n = 0:demilongZ
     k(i+1,j+1,n+1) = exp(-pas^2*(i^2+j^2+n^2)/lambda2);
    end
   end
  end
  for i = 2:demilongX
   for j = 2:demilongY
    for n = 2:demilongZ
     k(longX-i+2,longY-j+2,n) = k(i,j,n);
     k(longX-i+2,j,n) = k(i,j,n);
     k(i,longY-j+2,n) = k(i,j,n);
     k(i,j,longZ-n+2) = k(i,j,n);
     k(longX-i+2,longY-j+2,longZ-n+2) = k(i,j,n);
     k(longX-i+2,j,longZ-n+2) = k(i,j,n);
     k(i,longY-j+2,longZ-n+2) = k(i,j,n);
    end
   end
  end
  for i = 2:demilongX
   for j = 2:demilongY
    k(longX-i+2,j,1) = k(i,j,1);
    k(i,longY-j+2,1) = k(i,j,1);
    k(longX-i+2,longY-j+2,1) = k(i,j,1);
    k(longX-i+2,j,demilongZ+1) = k(i,j,demilongZ+1);
    k(i,longY-j+2,demilongZ+1) = k(i,j,demilongZ+1);
    k(longX-i+2,longY-j+2,demilongZ+1) = k(i,j,demilongZ+1);
   end
  end
  for i = 2:demilongX
   for n = 2:demilongZ
    k(longX-i+2,1,n) = k(i,1,n);
    k(i,1,longZ-n+2) = k(i,1,n);
    k(longX-i+2,1,longZ-n+2) = k(i,1,n);
    k(longX-i+2,demilongY+1,n) = k(i,demilongY+1,n);
    k(i,demilongY+1,longZ-n+2) = k(i,demilongY+1,n);
    k(longX-i+2,demilongY+1,longZ-n+2) = k(i,demilongY+1,n);
   end
  end
  for j = 2:demilongY
   for n = 2:demilongZ
    k(1,longY-j+2,n) = k(1,j,n);
    k(1,j,longZ-n+2) = k(1,j,n);
    k(1,longY-j+2,longZ-n+2) = k(1,j,n);
    k(demilongX+1,longY-j+2,n) = k(demilongX+1,j,n);
    k(demilongX+1,j,longZ-n+2) = k(demilongX+1,j,n);
    k(demilongX+1,longY-j+2,longZ-n+2) = k(demilongX+1,j,n);
   end
  end
  k((demilongX+2):longX,1,1) = k(demilongX:-1:2,1,1);
  k(1,(demilongY+2):longY,1) = k(1,demilongY:-1:2,1);
  k(1,1,(demilongZ+2):longZ) = k(1,1,demilongZ:-1:2);

  k((demilongX+2):longX,demilongY+1,1) = k(demilongX:-1:2,demilongY+1,1);
  k((demilongX+2):longX,1,demilongZ+1) = k(demilongX:-1:2,1,demilongZ+1);
  k((demilongX+2):longX,demilongY+1,demilongZ+1) = k(demilongX:-1:2,demilongY+1,demilongZ+1);

  k(1,(demilongY+2):longY,demilongZ+1) = k(1,demilongY:-1:2,demilongZ+1);
  k(demilongX+1,(demilongY+2):longY,1) = k(demilongX+1,demilongY:-1:2,1);
  k(demilongX+1,(demilongY+2):longY,demilongZ+1) = k(demilongX+1,demilongY:-1:2,demilongZ+1);

  k(demilongX+1,1,(demilongZ+2):longZ) = k(demilongX+1,1,demilongZ:-1:2);
  k(1,demilongY+1,(demilongZ+2):longZ) = k(1,demilongY+1,demilongZ:-1:2);
  k(demilongX+1,demilongY+1,(demilongZ+2):longZ) = k(demilongX+1,demilongY+1,demilongZ:-1:2);

  fft3k = fftn(k); % transformee de Fourier du noyau
  fft3k = real(fft3k); % pour éviter les imprécisions numériques
  fft3k_d = fft3k(1:(demilongX+1),:,:); % matrice decimée pour tenir compte des symmétries de fft3k
end
