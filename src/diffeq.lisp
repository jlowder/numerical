(in-package :cl-user)

(defpackage :numerical.diffeq
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+
                :genlambda/2)
  (:export :euler
           :modified-euler
           :improved-euler
           :runge-kutta/4
           :extrapolated-midpoint/4))

(in-package :numerical.diffeq)

(defmacro euler (fxy)
  "Create a lambda that can evaluate a first-order ordinary
differential equation `FXY`. The lambda takes the limits of the region
to evaluate as parameters, an `N` keyword argument for the number of
steps to take, and an `INITIAL` keyword argument for the initial
value. The lambda will use the Euler method to evaluate the result."
  (let ((gxy (gensym)))
    `(let ((,gxy (genlambda/2 ,fxy)))
       (lambda (a b &key (initial 0d0) (n 1))
         (let ((step (/ (- b a) n)))
           (loop repeat n
              as y = initial then yi
              as xi = a then (+ xi step)
              as m = (funcall ,gxy xi y)
              as yi = (+ y (* 1d0 step m))
              finally (return yi)))))))

(defmacro modified-euler (fxy)
  "Create a lambda that can evaluate a first-order ordinary
differential equation `FXY`. The lambda takes the limits of the region
to evaluate as parameters, an `N` keyword argument for the number of
steps to take, and an `INITIAL` keyword argument for the initial
value. The lambda will use the Modified Euler method (a.k.a. midpoint
method) to evaluate the result."
  (let ((gxy (gensym)))
    `(let ((,gxy (genlambda/2 ,fxy)))
       (lambda (a b &key (initial 0d0) (n 1))
         (let ((step (/ (- b a) n)))
           (loop repeat n
              as y = initial then yi
              as xi = a then (+ xi step)
              as m = (funcall ,gxy (+ xi (* 1/2 step)) y)
              as yi = (+ y (* 1d0 step m))
              finally (return yi)))))))

(defmacro improved-euler (fxy)
  "Create a lambda that can evaluate a first-order ordinary
differential equation `FXY`. The lambda takes the limits of the region
to evaluate as parameters, an `N` keyword argument for the number of
steps to take, and an `INITIAL` keyword argument for the initial
value. The lambda will use the Improved Euler method to evaluate the
result."
  (let ((gxy (gensym)))
    `(let ((,gxy (genlambda/2 ,fxy)))
       (lambda (a b &key (initial 0d0) (n 1))
         (let ((step (/ (- b a) n)))
           (loop repeat n
              as y = initial then yi
              as xi = a then (+ xi step)
              as f0 = (funcall ,gxy xi y)
              as m = (* 1/2 (+ f0 (funcall ,gxy (+ xi step) (+ y (* step f0)))))
              as yi = (+ y (* 1d0 step m))
              finally (return yi)))))))

(defmacro runge-kutta/4 (fxy)
  "Create a lambda that can evaluate a first-order ordinary
differential equation `FXY`. The lambda takes the limits of the region
to evaluate as parameters, an `N` keyword argument for the number of
steps to take, and an `INITIAL` keyword argument for the initial
value. The lambda will use the 4th Order Runge-Kutta method to
evaluate the result."
  (let ((gxy (gensym)))
    `(let ((,gxy (genlambda/2 ,fxy)))
       (lambda (a b &key (initial 0d0) (n 1))
         (let* ((step (/ (- b a) n))
                (h2 (/ step 2d0)))
           (loop repeat n
              as y = initial then yi
              as xi = a then (+ xi step)
              as f0 = (funcall ,gxy xi y)
              as f1 = (funcall ,gxy (+ xi h2) (+ y (* h2 f0)))
              as f2 = (funcall ,gxy (+ xi h2) (+ y (* h2 f1)))
              as f3 = (funcall ,gxy (+ xi h2) (+ y (* step f2)))
              as m = (* 1/6 (+ f0 f1 f1 f2 f2 f3))
              as yi = (+ y (* 1d0 step m))
              finally (return yi)))))))

(defmacro extrapolated-midpoint/4 (fxy)
  "Create a lambda that can evaluate a first-order ordinary
differential equation `FXY`. The lambda takes the limits of the region
to evaluate as parameters, an `N` keyword argument for the number of
steps to take, and an `INITIAL` keyword argument for the initial
value. The lambda will use the Modified Euler (i.e. midpoint) method
in combination with 4th order Richardson Extrapolation to evaluate the
result."
  (let ((gxy (gensym)))
    `(let ((,gxy (genlambda/2 ,fxy)))
       (labels ((extrapolate (l &optional (ni 2) nl)
                  (if (cadr l)
                      (let ((num (expt 2 ni)))
                        (extrapolate (cdr l) ni (append nl (list (/ (- (* num (cadr l)) (car l)) (1- num))))))
                      (if nl (extrapolate nl (+ ni 2)) (car l)))))
         (lambda (a b &key (initial 0d0) (n 1))
           (extrapolate (loop repeat 4
                           as m = 2 then (+ m m)
                           as step = (/ (- b a) (* m n))
                           collect (loop repeat (* n m)
                                      as y = initial then yi
                                      as xi = a then (+ xi step)
                                      as f0 = (funcall ,gxy xi y)
                                      as m = (* 1/2 (+ f0 (funcall ,gxy (+ xi step) (+ y (* step f0)))))
                                      as yi = (+ y (* 1d0 step m))
                                      finally (return yi)))))))))
