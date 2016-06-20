(in-package :cl-user)

(defpackage :numerical.test
  (:use :common-lisp 
        :lisp-unit 
        :numerical 
        :numerical.roots
        :numerical.interpolation
        :numerical.integration))

(in-package :numerical.test)

(defun pretty-close (a b)
  (< (abs (- a b)) 5d-7))

(defun really-close (a b)
  (< (abs (- a b)) 5d-15))

(defmacro ballpark (a b)
  `(if_ ,b
       (assert-true (< (abs (- ,a _)) 5d-3))
       (assert-equal ,a _)))

(defmacro ~= (a b)
  `(if_ ,b
        (assert-true (pretty-close ,a _))
        (assert-equal ,a _)))

(defmacro == (a b)
  `(if_ ,b
        (assert-true (really-close ,a _))
        (assert-equal ,a _)))

(defun !~= (a b)
  (assert-false (pretty-close a b)))

(defun !== (a b)
  (assert-false (really-close a b)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; root finding tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-test bisect
  (let ((exact (/ pi 2)))
    (~= exact (funcall (bisect #'cos) 0 3))
    (!== exact (funcall (bisect #'cos) 0 3))
    (== exact (funcall (bisect #'cos :tolerance 5d-15) 0 3))
    (~= exact (funcall (bisect+ (x)
                                (cos x))
                       0 3))
    (!== exact (funcall (bisect+ (x)
                                 (cos x))
                        0 3))
    (== exact (funcall (bisect+ (x :tolerance 5d-15)
                                (cos x))
                       0 3))))

;; all the other root-finders produce "exact" answers without having to raise the tolerance
(define-test false-position
  (let ((exact (/ pi 2)))
    (== exact (funcall (false-position #'cos) 0 3))
    (== exact (funcall (false-position+ (x)
                                        (cos x))
                       0 3))))

(define-test secant/bisect
  (let ((exact (/ pi 2)))
    (== exact (funcall (secant/bisect #'cos) 0 3))
    (== exact (funcall (secant/bisect+ (x)
                                       (cos x))
                       0 3))))

(define-test quadratic/bisect
  (let ((exact (/ pi 2)))
    (== exact (funcall (quadratic/bisect #'cos) 0 3))
    (== exact (funcall (quadratic/bisect+ (x)
                                          (cos x))
                       0 3))))

(define-test newton-raphson
  (let ((exact (/ pi 2)))
    (== exact (funcall (newton-raphson #'cos (lambda (x) (* -1 (sin x)))) .5))
    (== exact (funcall (newton-raphson+ (x)
                                        (cos x)
                                        (* -1 (sin x)))
                       1/2))
    (== exact (funcall (newton-raphson/bisect #'cos (lambda (x) (* -1 (sin x)))) 0 3))
    (== exact (funcall (newton-raphson/bisect+ (x)
                                               (cos x)
                                               (* -1 (sin x)))
                       0 3))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; interpolation tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun j (n)
  (nth n '((1.0d0 0.7651976866d0 0.2238907791d0 -0.2600519549d0
            -0.3971498099d0 -0.1775967713d0 0.1506452573d0
            0.3000792705d0 0.1716508071d0 -0.0903336112d0
            -0.2459357645d0)
           (0.0d0 0.4400505857d0 0.5767248078d0 0.3390589585d0
            -0.0660433280d0 -0.3275791376d0 -0.2766838581d0
            -0.0046828235d0 0.2346363469d0 0.2453117866d0
            0.0434727462d0)
           (0d0 0.1149034849d0 0.3528340286d0 0.4860912606d0
            0.3641281459d0 0.0465651163d0 -0.2428732100d0
            -0.3014172201d0 -0.1129917204d0 0.1448473415d0
            0.2546303137d0))))

(defparameter j1 (j 1))

(define-test lagrange-polynomial
  (let ((interpolator (lagrange-polynomial (loop for i from 0 to 10 collect i) j1)))
    (== 0.5579736040706821d0 (funcall interpolator 1.5))
    (== 0.5559568307446474d0 (funcall interpolator 2.2))
    (== 0.06790220312575382d0 (funcall interpolator 9.9))
    (== 0.24194621903592411d0 (funcall interpolator .5))))

(define-test lagrange/order
  (let ((interpolator (lagrange/order (loop for i from 0 to 10 collect i) j1 5)))
    (== 0.557257489765625d0 (funcall interpolator 1.5))
    (== 0.556519219757108d0 (funcall interpolator 2.2))
    (== 0.06837787882385432d0 (funcall interpolator 9.9))
    (== 0.24265791351562496d0 (funcall interpolator .5))))

(define-test cubic-spline
  (let ((interpolator (cubic-spline (loop for i from 0 to 10 collect i) j1)))
    (== 0.5568778304350899d0 (funcall interpolator 1.5))
    (== 0.5557082105761957d0 (funcall interpolator 2.2))
    (== 0.06780741900679668d0 (funcall interpolator 9.9))
    (== 0.24178396040497008d0 (funcall interpolator .5))))

(define-test hermite/cubic
  (let ((interpolator (hermite/cubic (loop for i from 0 to 10 collect i)
                                     j1
                                     (loop for i from 0 to 10 collect (* 1/2 (- (nth i (j 0)) (nth i (j 2))))))))
    (== 0.5570900374500001d0 (funcall interpolator 1.5))
    (== 0.5556934847566745d0 (funcall interpolator 2.2))
    (== 0.06833894614750796d0 (funcall interpolator 9.9))
    (== 0.24188190524375d0 (funcall interpolator .5))))

(define-test cubic-spline-prime
  (let ((differentiator (cubic-spline-prime (loop for i from 0 to 10 collect i) j1))
        (exact 1.8484212499354509d0))
    (== exact (funcall (quadratic/bisect+ (x)
                                          (funcall differentiator x)) 1 3))))

(define-test first-derivatives
  (labels ((f (x)
             (* x (exp x)))
           (df (x)
             (+ (f x) (exp x)))
           (ddf (x)
             (+ (f x) (exp x) (exp x))))
    (let ((exact (df 2d0))
          (d1 (first-derivative/2 #'f))
          (d2 (first-derivative/4 #'f))
          (d3 (first-derivative #'f :tolerance 5d-15))
          (d3+ (first-derivative+ (x :tolerance 5d-15)
                                  (* x (exp x)))))
      (!~= exact (funcall d1 2))
      (ballpark exact (funcall d2 2))
      (~= exact (funcall d3 2))
      (~= exact (funcall d3+ 2))
      )))

(define-test second-derivatives
  (labels ((f (x)
             (* x (exp x)))
           (df (x)
             (+ (f x) (exp x)))
           (ddf (x)
             (+ (f x) (exp x) (exp x))))
    (let ((exact (ddf 2d0))
          (d1 (second-derivative/3 #'f))
          (d2 (second-derivative/5 #'f))
          (d3 (second-derivative #'f :step .4 :tolerance 5d-15))
          (d3+ (second-derivative+ (x :step .4 :tolerance 5d-15)
                                   (* x (exp x)))))
      (!~= exact (funcall d1 2))
      (ballpark exact (funcall d2 2))
      (ballpark exact (funcall d3 2)) ;; todo figure out why not more precise
      (ballpark exact (funcall d3+ 2)) ;; todo figure out why not more precise
      )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; integration tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-test piecewise/2
  ;; piecewise/2 should be exact for a straight line
  (labels ((f (x) (* 2 x))
           (exact (a b) (- (* b b) (* a a))))
    (let ((exact (exact 1 5))
          (p2 (piecewise/2 #'f))
          (p2+ (piecewise/2+ (x) (* 2 x))))
      (== exact (funcall p2 1 5))
      (== exact (funcall p2+ 1 5)))))

(define-test piecewise/3
  ;; piecewise/3 should be exact for a quadratic line
  (labels ((f (x) (* 3 x x))
           (exact (a b) (- (* b b b) (* a a a))))
    (let ((exact (exact 1 5))
          (p3 (piecewise/3 #'f))
          (p3+ (piecewise/3+ (x) (* 3 x x))))
      (== exact (funcall p3 1 5))
      (== exact (funcall p3+ 1 5)))))

(define-test piecewise/4
  ;; piecewise/4 should be exact for a cubic line
  (labels ((f (x) (* 4 x x x))
           (exact (a b) (- (* b b b b) (* a a a a))))
    (let ((exact (exact 1 5))
          (p4 (piecewise/4 #'f))
          (p4+ (piecewise/4+ (x) (* 4 x x x))))
      (== exact (funcall p4 1 5))
      (== exact (funcall p4+ 1 5)))))

(define-test piecewise/5
  ;; piecewise/5 should be exact for a quartic line
  (labels ((f (x) (* 5 x x x x))
           (exact (a b) (- (* b b b b b) (* a a a a a))))
    (let ((exact (exact 1 5))
          (p5 (piecewise/5 #'f))
          (p5+ (piecewise/5+ (x) (* 5 x x x x))))
      (== exact (funcall p5 1 5))
      (== exact (funcall p5+ 1 5)))))

(define-test romberg
  (labels ((f (x) (* 5 x x x x))
           (exact (a b) (- (* b b b b b) (* a a a a a))))
    (let ((exact (exact 1 5))
          (r (romberg #'f))
          (r+ (romberg+ (x)
                        (* 5 x x x x))))
      (== exact (funcall r 1 5))
      (== exact (funcall r+ 1 5)))))

(define-test euler-mcclaurin
  (labels ((f (x) (* 5 x x x x))
           (df (x) (* 20 x x x))
           (df2 (x) (* 60 x x))
           (df3 (x) (* 120 x))
           (exact (a b) (- (* b b b b b) (* a a a a a))))
    (let ((exact (exact 1 5))
          (em (euler-mcclaurin #'f #'df #'df3))
          (em+ (euler-mcclaurin+ (x)
                                 (* 5 x x x x)
                                 (* 20 x x x)
                                 (* 120 x))))
      (== exact (funcall em 1 5))
      (== exact (funcall em+ 1 5)))))

(define-test gauss-legendre/2
  (labels ((f (x) (* 3 x x))
           (exact (a b) (- (* b b b) (* a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/2 #'f))
          (gl+ (gauss-legendre/2+ (x)
                                   (* 3 x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/3
  (labels ((f (x) (* 4 x x x))
           (exact (a b) (- (* b b b b) (* a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/3 #'f))
          (gl+ (gauss-legendre/3+ (x)
                                  (* 4 x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/4
  (labels ((f (x) (* 4 x x x))
           (exact (a b) (- (* b b b b) (* a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/4 #'f))
          (gl+ (gauss-legendre/4+ (x)
                                  (* 4 x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/5
  (labels ((f (x) (* 5 x x x x))
           (exact (a b) (- (* b b b b b) (* a a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/5 #'f))
          (gl+ (gauss-legendre/5+ (x)
                                  (* 5 x x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/6
  (labels ((f (x) (* 6 x x x x x))
           (exact (a b) (- (* b b b b b b) (* a a a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/6 #'f))
          (gl+ (gauss-legendre/6+ (x)
                                  (* 6 x x x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/7
  (labels ((f (x) (* 7 x x x x x x))
           (exact (a b) (- (* b b b b b b b) (* a a a a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/7 #'f))
          (gl+ (gauss-legendre/7+ (x)
                                  (* 7 x x x x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/8
  (labels ((f (x) (* 8 x x x x x x x))
           (exact (a b) (- (* b b b b b b b b) (* a a a a a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/8 #'f))
          (gl+ (gauss-legendre/8+ (x)
                                  (* 8 x x x x x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(define-test gauss-legendre/9
  (labels ((f (x) (* 9 x x x x x x x x))
           (exact (a b) (- (* b b b b b b b b b) (* a a a a a a a a a))))
    (let ((exact (exact 1 5))
          (gl (gauss-legendre/9 #'f))
          (gl+ (gauss-legendre/9+ (x)
                                  (* 9 x x x x x x x x))))
      (~= exact (funcall gl 1 5))
      (~= exact (funcall gl+ 1 5)))))

(let ((*print-errors* t)
      (*print-failures* t))
  (run-tests))
