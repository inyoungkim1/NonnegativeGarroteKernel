strassenInv <- function(A){
   if (nrow(A) != ncol(A))
       { stop("only square matrices can be inverted") }
   is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)
   abs(x - round(x)) < tol

   if ( (is.wholenumber(log(nrow(A), 2)) != TRUE) || (is.wholenumber(log(ncol(A), 2)) != TRUE) )
       { stop("only square matrices of dimension 2^k * 2^k can be inverted with Strassen method") }

 A <- div4(A, dim(A)[1])
 R1 <- solve(A$X11)
 R2 <- A$X21 %*% R1
 R3 <- R1 %*% A$X12
 R4 <- A$X21 %*% R3
 R5 <- R4 - A$X22
 R6 <- solve(R5)
 C12 <- R3 %*% R6
 C21 <- R6 %*% R2
 R7 <- R3 %*% C21
 C11 <- R1 - R7
 C22 <- -R6
 C <- rbind(cbind(C11,C12), cbind(C21,C22))
 return(C)}
 
 
strassenInv2 <- function(A){
   if (nrow(A) != ncol(A))
     { stop("only square matrices can be inverted") }
   is.wholenumber <-
     function(x, tol = .Machine$double.eps^0.5)
   abs(x - round(x)) < tol
   if ( (is.wholenumber(log(nrow(A), 2)) != TRUE) || (is.wholenumber(log(ncol(A), 2)) != TRUE) )
     { stop("only square matrices of dimension 2^k * 2^k can be inverted with Strassen method") }

   A <- div4(A, dim(A)[1])
   R1 <- strassenInv(A$X11)
   R2 <- strassen(A$X21 , R1)
   R3 <- strassen(R1 , A$X12)
   R4 <- strassen(A$X21 , R3)
   R5 <- R4 - A$X22
   R6 <- strassenInv(R5)
   C12 <- strassen(R3 , R6)
   C21 <- strassen(R6 , R2)
   R7 <- strassen(R3 , C21)
   C11 <- R1 - R7
   C22 <- -R6
   C <- rbind(cbind(C11,C12), cbind(C21,C22))

return(C)}


strassenInv3 <- function(A){
    if (nrow(A) != ncol(A))
       { stop("only square matrices can be inverted") }
    is.wholenumber <-
       function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if ( (is.wholenumber(log(nrow(A), 2)) != TRUE) || (is.wholenumber(log(ncol(A), 2)) != TRUE) )
       { stop("only square matrices of dimension 2^k * 2^k can be inverted with Strassen method") }

 A <- div4(A, dim(A)[1])
 R1 <- strassenInv2(A$X11)
 R2 <- strassen2(A$X21 , R1)
 R3 <- strassen2(R1 , A$X12)
 R4 <- strassen2(A$X21 , R3)
 R5 <- R4 - A$X22
 R6 <- strassenInv2(R5)
 C12 <- strassen2(R3 , R6)
 C21 <- strassen2(R6 , R2)
 R7 <- strassen2(R3 , C21)
 C11 <- R1 - R7
 C22 <- -R6

 C <- rbind(cbind(C11,C12), cbind(C21,C22))
 return(C)
 }

################# multiplicarion ###################
div4 <- function(A, r){
  A <- list(A)
  A11 <- A[[1]][1:(r/2),1:(r/2)]
  A12 <- A[[1]][1:(r/2),(r/2+1):r]
  A21 <- A[[1]][(r/2+1):r,1:(r/2)]
  A22 <- A[[1]][(r/2+1):r,(r/2+1):r]
  A <- list(X11=A11, X12=A12, X21=A21, X22=A22)
  return(A)
 }

 strassen <- function(A, B){
  A <- div4(A, dim(A)[1])
  B <- div4(B, dim(B)[1])
  M1 <- (A$X11+A$X22) %*% (B$X11+B$X22)
  M2 <- (A$X21+A$X22) %*% B$X11
  M3 <- A$X11 %*% (B$X12-B$X22)
  M4 <- A$X22 %*% (B$X21-B$X11)
  M5 <- (A$X11+A$X12) %*% B$X22
  M6 <- (A$X21-A$X11) %*% (B$X11+B$X12)
  M7 <- (A$X12-A$X22) %*% (B$X21+B$X22)

  C11 <- M1+M4-M5+M7
  C12 <- M3+M5
  C21 <- M2+M4
  C22 <- M1-M2+M3+M6

  C <- rbind(cbind(C11,C12), cbind(C21,C22))
  return(C)
 }



 strassen2 <- function(A, B){
  A <- div4(A, dim(A)[1])
  B <- div4(B, dim(B)[1])
  M1 <- strassen((A$X11+A$X22) , (B$X11+B$X22))
  M2 <- strassen((A$X21+A$X22) , B$X11)
  M3 <- strassen(A$X11 , (B$X12-B$X22))
  M4 <- strassen(A$X22 , (B$X21-B$X11))
  M5 <- strassen((A$X11+A$X12) , B$X22)
  M6 <- strassen((A$X21-A$X11) , (B$X11+B$X12))
  M7 <- strassen((A$X12-A$X22) , (B$X21+B$X22))

  C11 <- M1+M4-M5+M7
  C12 <- M3+M5
  C21 <- M2+M4
  C22 <- M1-M2+M3+M6

  C <- rbind(cbind(C11,C12), cbind(C21,C22))
  return(C)
 }
 






