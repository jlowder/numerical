(in-package :cl-user)

(defpackage :numerical.diffeq
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+)
  (:export :euler
           :modified-euler
           :improved-euler
           :runge-kutta/4))

(in-package :numerical.diffeq)

(defun euler (fxy)
  (lambda (a b &key (initial 0d0) (n 1))
    (let ((step (/ (- b a) n)))
      (loop repeat n
         as y = initial then yi
         as xi = a then (+ xi step)
         as m = (funcall fxy xi y)
         as yi = (+ y (* 1d0 step m))
         finally (return yi)))))

(defun modified-euler (fxy)
  (lambda (a b &key (initial 0d0) (n 1))
    (let ((step (/ (- b a) n)))
      (loop repeat n
         as y = initial then yi
         as xi = a then (+ xi step)
         as m = (funcall fxy (+ xi (* 1/2 step)) y)
         as yi = (+ y (* 1d0 step m))
         finally (return yi)))))

(defun improved-euler (fxy)
  (lambda (a b &key (initial 0d0) (n 1))
    (let ((step (/ (- b a) n)))
      (loop repeat n
         as y = initial then yi
         as xi = a then (+ xi step)
         as f0 = (funcall fxy xi y)
         as m = (* 1/2 (+ f0 (funcall fxy (+ xi step) (+ y (* step f0)))))
         as yi = (+ y (* 1d0 step m))
         finally (return yi)))))

(defun runge-kutta/4 (fxy)
  (lambda (a b &key (initial 0d0) (n 1))
    (let* ((step (/ (- b a) n))
           (h2 (/ step 2d0)))
      (loop repeat n
         as y = initial then yi
         as xi = a then (+ xi step)
         as f0 = (funcall fxy xi y)
         as f1 = (funcall fxy (+ xi h2) (+ y (* h2 f0)))
         as f2 = (funcall fxy (+ xi h2) (+ y (* h2 f1)))
         as f3 = (funcall fxy (+ xi h2) (+ y (* step f2)))
         as m = (* 1/6 (+ f0 f1 f1 f2 f2 f3))
         as yi = (+ y (* 1d0 step m))
         finally (return yi)))))

