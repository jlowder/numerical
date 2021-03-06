(in-package :cl-user)

(defpackage :numerical.roots
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+
                :genlambda/1)
  (:export :bisect
           :newton-raphson
           :newton-raphson/bisect
           :root-brackets
           :false-position
           :secant/bisect
           :quadratic/bisect))

(in-package :numerical.roots)

(defmacro bisect (f &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for finding a root of `F` between two values, to a
specified accuracy and within a specified number of iteration
attempts. The lambda takes the lower and upper values as parameters
and will use a brute-force bisection algorithm to find a root."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (labels ((rec (left right &optional (cnt 0) (fleft (funcall ,g left)) (fright (funcall ,g right)))
                  (let* ((middle (* 0.5d0 (+ left right)))
                         (fmid (funcall ,g middle)))
                    (if (or
                         (> (* fleft fright) 0d0) ; not bracketed
                         (>= cnt ,limit))
                        (values middle nil cnt)
                        (if (< (abs (/_ (- right left) middle)) ,tolerance)
                            (values middle t cnt)
                            (if (< (* fleft fmid) 0)
                                (rec left middle (1+ cnt) fleft fmid)
                                (rec middle right (1+ cnt) fmid fright)))))))
         #'rec))))

(defmacro newton-raphson (f df &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for finding a root of `F`, to a specified accuracy
and within a specified number of iteration attempts. The lambda will
use the newton-raphson algorithm (which requires the derivative of the
function to be provided as `DF`) to converge on the solution. The
lambda takes a guess to the solution as a parameter."
  (let ((g (gensym))
        (dg (gensym)))
    `(let ((,g (genlambda/1 ,f))
           (,dg (genlambda/1 ,df)))
       (labels ((rec (guess &optional (cnt 0))
                  (if (>= cnt ,limit)
                      (values guess nil ,limit)
                      (let ((delta (/_ (funcall ,g guess) (funcall ,dg guess) -1d0)))
                        (if (<= (abs (/_ delta (+ delta guess))) ,tolerance)
                            (values (+ guess delta) t cnt)
                            (rec (+ guess delta) (1+ cnt)))))))
         #'rec))))

(defmacro newton-raphson/bisect (f df &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for finding a root of `F` between two values, to a
specified accuracy and within a specified number of iteration
attempts. The lambda will use a hybrid algorithm that decides on each
iteration whether to take a newton-raphson step or a bisection
step. The lambda takes the lower and upper values as parameters."
  (let ((g (gensym))
        (dg (gensym)))
    `(let ((,g (genlambda/1 ,f))
           (,dg (genlambda/1 ,df)))
       (labels ((rec (left right &optional
                           (cnt 0)
                           (fleft (funcall ,g left))
                           (fright (funcall ,g right))
                           (best (if (<= (abs fleft) (abs fright)) left right))
                           (fbest (funcall ,g best))
                           (dfbest (funcall ,dg best))
                           (nr (lambda () (let ((delta (/_ fbest dfbest -1d0)))
                                            (values delta (+ best delta)))))
                           (bisect (lambda () (values (/_ (- right left) 2d0)
                                                      (/_ (+ right left) 2d0)))))
                  (if (or
                       (> (* fleft fright) 0d0) ; not bracketed
                       (> cnt ,limit)) ; fail to converge
                      (values best nil cnt)
                      (multiple-value-bind (delta best) (funcall (if (<= (*
                                                                          (- (* dfbest (- best left)) fbest)
                                                                          (- (* dfbest (- best right)) fbest))
                                                                         0d0)
                                                                     nr bisect))
                        (if (< (abs (/_ delta best)) ,tolerance)
                            (values best t cnt)
                            (if (>= (* (funcall ,g best) fleft) 0d0)
                                (rec best right (1+ cnt))
                                (rec left best (1+ cnt))))))))
         #'rec))))

(defmacro root-brackets (f &key (steps 100))
  "Create a lambda that will search for pairs of values that bracket
the roots of `F`. `STEPS` determines the number of intervals that will
be sampled. A list of value pairs will be returned."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (lambda (left right)
         (loop
            for px = left then x
            for x = left then (+ x (/ (- right left) ,steps))
            when (not (equal (plusp (funcall ,g x))
                             (plusp (funcall ,g px))))
            collect (list (coerce px 'double-float) (coerce x 'double-float)) into res
            until (>= x right)
            finally (return res))))))

(defmacro false-position (f &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for attempting to find a root of `F` between two
values, to a specified accuracy and within a specified number of
iteration attempts. The lambda will use the method of false position
to converge on the solution. The lambda takes the lower and upper
values as parameters."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (labels ((rec (left right &optional
                           (cnt 0)
                           (prev right)
                           (fleft (funcall ,g left))
                           (fright (funcall ,g right)))
                  (if (> cnt ,limit)
                      (values prev nil cnt)
                      (let* ((guess (/_ (- (* left fright) (* 1d0 right fleft))
                                        (- fright fleft)))
                             (fguess (funcall ,g guess))
                             (err (abs (/_ (- guess prev) guess))))
                        (if (<= err ,tolerance)
                            (values guess t cnt)
                            (if (< (* fleft fguess) 0d0)
                                (rec left guess (1+ cnt) guess fleft fguess)
                                (rec guess right (1+ cnt) guess fguess fright)))))))
         #'rec))))

(defmacro secant/bisect (f &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for attempting to find a root of `F` between two
values, to a specified accuracy and within a specified number of
iteration attempts. The lambda will use the Secant method to converge
on the solution, but reverts to the bisection method if unable to
converge. The lambda takes the lower and upper values as parameters."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (let ((bistepper (bisect ,f :tolerance ,tolerance :limit ,1)))
         (labels ((secant (xi-1 xi)
                    (let ((fxi-1 (funcall ,g xi-1))
                          (fxi (funcall ,g xi)))
                      (- xi (* 1d0 fxi (/ (- xi xi-1)
                                          (- fxi fxi-1))))))
                  (rec (left right &optional
                             (cnt 0)
                             (fleft (funcall ,g left)) 
                             (fright (funcall ,g right))
                             (best (if (<= (abs fleft) (abs fright)) left right))
                             (xi-1 (funcall bistepper left right 1))
                             (xi (funcall bistepper left right 0)))
                    (if (> (* fleft fright) 0d0)
                        (values best nil cnt)
                        (if (> cnt ,limit)
                            (funcall (bisect ,f :tolerance ,tolerance :limit ,limit) left right)
                            (let* ((guess (secant xi-1 xi))
                                   (fguess (funcall ,g guess))
                                   (err (abs (/_ (- guess xi) guess))))
                              (if (<= err ,tolerance)
                                  (values guess t cnt)
                                  (if (<= (* fleft fguess) 0d0)
                                      (rec left guess (1+ cnt) fleft fguess guess xi guess)
                                      (rec guess right (1+ cnt) fguess fright guess xi guess))))))))
           #'rec)))))

(defmacro quadratic/bisect (f &key (tolerance +tolerance+) (limit +limit+))
  "Create a lambda for attempting to find a root of `F` between two
values, to a specified accuracy and within a specified number of
iteration attempts. The lambda will use a quadratic approximation
method to converge on the solution, but reverts to the bisection
method if unable to converge. The lambda takes the lower and upper
values as parameters."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (let ((bistepper (bisect ,f :tolerance ,tolerance :limit 2)))
         (labels ((quadratic (x0 x1 x2)
                    (if (or
                         (equal x0 x1)
                         (equal x1 x2)
                         (equal x0 x2))
                        x2
                        (let* ((fx0 (funcall ,g x0))
                               (fx1 (funcall ,g x1))
                               (fx2 (funcall ,g x2))
                               (c fx2)
                               (b (/_ (- (* (- x0 x2) (- x0 x2) (- fx1 fx2)) (* (- x1 x2) (- x1 x2) (- fx0 fx2)))
                                      (* (- x0 x1) (- x0 x2) (- x1 x2))))
                               (a (/_ (- (* (- x1 x2) (- fx0 fx2)) (* (- x0 x2) (- fx1 fx2)))
                                      (* (- x0 x1) (- x0 x2) (- x1 x2))))
                               (numer (- 0d0 c c))
                               (denom+ (+ b (expt (- (* b b) (* 4 a c)) .5d0)))
                               (denom- (- b (expt (- (* b b) (* 4 a c)) .5d0))))
                          (+ x2 (if (>= (abs denom+) (abs denom-))
                                    (/_ numer denom+)
                                    (/_ numer denom-))))))
                  (rec (left right &optional
                             (cnt 0)
                             (fleft (funcall ,g left)) 
                             (fright (funcall ,g right))
                             (best (if (<= (abs fleft) (abs fright)) left right))
                             (xi-2 (funcall bistepper left right 2))
                             (xi-1 (funcall bistepper left right 1))
                             (xi (funcall bistepper left right 0)))
                    (if (> (* fleft fright) 0d0)
                        (values best nil cnt)
                        (if (> cnt ,limit)
                            (funcall (bisect ,f :tolerance ,tolerance :limit ,limit) left right)
                            (let* ((guess (quadratic xi-2 xi-1 xi))
                                   (fguess (funcall ,g guess))
                                   (err (abs (/_ (- guess xi) guess))))
                              (if (<= err ,tolerance)
                                  (values guess t cnt)
                                  (if (<= (* fleft fguess) 0d0)
                                      (rec left guess (1+ cnt) fleft fguess guess xi-1 xi guess)
                                      (rec guess right (1+ cnt) fguess fright guess xi-1 xi guess))))))))
           #'rec)))))
