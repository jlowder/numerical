(in-package :cl-user)

(defpackage :numerical.integration
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+
                :genlambda/1)
  (:export :piecewise/2
           :piecewise/3
           :piecewise/4
           :piecewise/5
           :romberg
           :euler-mcclaurin
           :gauss-laguerre/2
           :gauss-laguerre/4
           :gauss-laguerre/6
           :gauss-laguerre/8
           :gauss-legendre/2
           :gauss-legendre/3
           :gauss-legendre/4
           :gauss-legendre/5
           :gauss-legendre/6
           :gauss-legendre/7
           :gauss-legendre/8
           :gauss-legendre/9
           :gauss-legendre/n
           :gauss-legendre
           :monte-carlo
           :standard-deviation
           :running-std-dev
           :running-monte-carlo))

(in-package :numerical.integration)

(defmacro piecewise/2 (f)
  "Create a lambda that integrates a region over `F`. The lambda takes
the lower and upper limits of the region as parameters, and optionally
the number of segments to divide the region into. The algorithm used
is piecewise linear approximation, i.e. the trapezoid rule."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (lambda (l u &optional (n 1))
         (let ((h (/ (- u l) n)))
           (* .5d0 h (loop repeat (1+ n)
                        as x = l then (+ x h)
                        as i = 0 then (1+ i)
                        as fx = (funcall ,g x)
                        summing (cond ((or (eq i 0) (eq i n)) fx)
                                      (t (* 2 fx))))))))))

(defmacro piecewise/3 (f)
  "Create a lambda that integrates a region over `F`. The lambda takes
the lower and upper limits of the region as parameters, and optionally
the number of segments to divide the region into (which will be
incremented to make even if necessary). The algorithm used is
piecewise quadratic approximation, i.e. Simpson's 1/3 rule."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (lambda (l u &optional (n 2))
         (let* ((n (if (oddp n) (1+ n) n)) ; force n to be even
                (h (/ (- u l) n)))
           (* (/ 1d0 3) h (loop repeat (1+ n)
                             as x = l then (+ x h)
                             as i = 0 then (1+ i)
                             as fx = (funcall ,g x)
                             summing (cond ((or (eq i 0) (eq i n)) fx)
                                           ((oddp i) (* 4 fx))
                                           (t (* 2 fx))))))))))

(defmacro piecewise/4 (f)
  "Create a lambda that integrates a region over `F`. The lambda takes
the lower and upper limits of the region as parameters, and optionally
the number of segments to divide the region into (which will be
incremented to make divisible by 3 if necessary). The algorithm used is
piecewise cubic approximation, i.e. Simpson's 3/8 rule."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (lambda (l u &optional (n 3))
         (let* ((n (- n (second (multiple-value-list (ceiling n 3))))) ; multiple of 3
                (h (/ (- u l) n)))
           (* .375d0 h (loop repeat (1+ n)
                          as x = l then (+ x h)
                          as i = 0 then (1+ i)
                          as fx = (funcall ,g x)
                          summing (cond ((or (eq i 0) (eq i n)) fx)
                                        ((eq (mod i 3) 0) (* 2 fx))
                                        (t (* 3 fx))))))))))

(defmacro piecewise/5 (f)
  "Create a lambda that integrates a region over `F`. The lambda takes
the lower and upper limits of the region as parameters, and optionally
the number of segments to divide the region into (which will be
incremented to make divisible by 4 if necessary). The algorithm used is
piecewise quartic approximation, i.e. Boole's (or Bode's) rule."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (lambda (l u &optional (n 4))
         (let* ((n (- n (second (multiple-value-list (ceiling n 4))))) ; multiple of 4
                (h (/ (- u l) n)))
           (* (/ 2d0 45) h (loop repeat (1+ n)
                              as x = l then (+ x h)
                              as i = 0 then (1+ i)
                              as fx = (funcall ,g x)
                              summing (cond ((or (eq i 0) (eq i n)) (* 7 fx))
                                            ((eq (mod i 4) 0) (* 14 fx))
                                            ((eq (mod i 4) 2) (* 12 fx))
                                            (t (* 32 fx))))))))))

(defmacro romberg (f &key (integrator 'piecewise/2) ratio-limit (tolerance +tolerance+) (limit +limit+))
  "Create a lambda that attempts to integrate a region over `F` to the
specified accuracy, within the specified limit of iteration, by
applying richardson extrapolation to the specified integrator. If
`RATIO-LIMIT` is specified, it is the amount of deviation allowed to
the ratio of error reduction between consecutive steps, which should
normally be around 4; the deviation allowed should be around 0.4, or
leave NIL to disable this feature."
  `(let ((integrator (,integrator ,f)))
     (labels ((rec (l u &optional (i 1) (Ih*2 (funcall integrator l u 1)) Pext Ih*4)
                (let* ((n (expt 2 i))
                       (Ih (funcall integrator l u n))
                       (div (- (expt 2d0 (1+ i)) 1))
                       (err3 (* 1/3 (- Ih Ih*2)))
                       (ext (+ Ih err3))
                       (err (if Pext (* (/_ 1d0 div) (- ext Pext)) err3))
                       (ext2 (+ ext err))
                       (R (if (and ,ratio-limit Ih*4 (not (equal Ih*2 Ih))) (abs (/_ (- Ih*4 Ih*2)
                                                                                     (- Ih*2 Ih)))
                              4)))
                  (if (and Pext
                           (> i 3)
                           (not (zerop ext2))
                           (< (abs (/_ (- Ih Ih*2) Ih)) ,tolerance))
                      (values ext2 t i)
                      (if (and (< i ,limit)
                               (or (<= i 3) (> (if ,ratio-limit ,ratio-limit .4d0) (abs (- 4 R)))))
                          (rec l u (1+ i) Ih ext2 Ih*2)
                          (values ext2 nil i))))))
       #'rec)))

(defmacro euler-mcclaurin (f &optional (df nil df-p) (df3 nil df3-p))
  "Create a lambda that integrates a region over `F`, where `DF` is
the first derivative function of `F` and `DF3` is the third derivative
function of `F`. `DF` and `DF3` are both optional; if neither is
specified, the behavior is identical to piecewise/2. Otherwise the
derivative functions are used with the Euler-McClaurin integration
rule to improve accuracy. The lambda takes the lower and upper limits
of the region as parameters, and optionally the number of segments to
divide the region into."
  (let ((dg (gensym))
        (dg3 (gensym)))
    `(let ((,dg (if ,df-p (genlambda/1 ,df) nil))
           (,dg3 (if ,df3-p (genlambda/1 ,df3) nil)))
       (lambda (l u &optional (n 1) (integrator (piecewise/2 ,f)))
         (let ((h (/ (- u l) n))
               (df (if ,df-p
                       (- (funcall ,dg l) (funcall ,dg u))
                       0))
               (df3 (if ,df3-p
                        (- (funcall ,dg3 l) (funcall ,dg3 u))
                        0)))
           (+ (funcall integrator l u n) (* (/ 1d0 12) h h df) (* (/ -1d0 720) h h h h df3)))))))

(defmacro gaussian-quadrature/n (f n)
  "Create a lambda that integrates a region over `F`. The lambda takes
the lower and upper limits of the region as parameters. `N` specifies
the order of the Legendre polynomial to use for the integration rule,
which can range from 2 to 9."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (flet ((gq (ti ci)
                (lambda (a b)
                  (let ((m (/ (- b a) 2d0))
                        (c (/ (+ b a) 2d0)))
                    (flet ((Ft (x)
                             (funcall ,g (+ c (* m x)))))
                      (* m (loop for v in ti
                              for z in ci
                              summing (* z (Ft v)))))))))
         (apply #'gq (cond ((eql ,n 2) '((-.577350269189626d0 .577350269189626d0)
                                         (1 1)))
                           ((eql ,n 3) '((-0.774596669241483d0 0d0 0.774596669241483d0)
                                         (5/9 8/9 5/9)))
                           ((eql ,n 4) '((-0.861136311594053d0 -0.339981043584856d0 0.339981043584856d0 0.861136311594053d0)
                                         (0.347854845137454d0 0.652145154862546d0 0.652145154862546d0 0.347854845137454d0)))
                           ((eql ,n 5) '((-0.906179845938664d0 -0.538469310105683d0 0d0 0.538469310105683d0 0.906179845938664d0)
                                         (0.236926885056189d0 0.478628670499366d0 0.5688888888889d0 0.478628670499366d0 0.236926885056189d0)))
                           ((eql ,n 6) '((-0.932469514203152d0 -0.661209386466265d0 -0.238619186083197d0
                                          0.238619186083197d0 0.661209386466265d0 0.932469514203152d0)
                                         (0.171324492379170d0 0.360761573048139d0 0.467913934572691d0
                                          0.467913934572691d0 0.360761573048139d0 0.171324492379170d0)))
                           ((eql ,n 7) '((-0.949107912342759d0 -0.741531185599394d0 -0.405845151377397d0 0d0 0.405845151377397d0
                                          0.741531185599394d0 0.949107912342759d0)
                                         (0.129484966168870d0 0.279705391489277d0 0.381830050505119d0 0.417959183673469d0
                                          0.381830050505119d0 0.279705391489277d0 0.129484966168870d0)))
                           ((eql ,n 8) '((-0.960289856497536d0 -0.796666477413627d0 -0.525532409916329d0 -0.183434642495650d0
                                          0.183434642495650d0 0.525532409916329d0 0.796666477413627d0 0.960289856497536d0)
                                         (0.101228536290376d0 0.222381034453374d0 0.313706645877887d0 0.362683783378362d0
                                          0.362683783378362d0 0.313706645877887d0 0.222381034453374d0 0.101228536290376d0)))
                           ((eql ,n 9) '((-0.968160239507626d0 -0.836031107326636d0 -0.613371432700590d0 -0.324253423403809d0 0d0
                                          0.324253423403809d0 0.613371432700590d0 0.836031107326636d0 0.968160239507626d0)
                                         (0.081274388361574d0 0.180648160694857d0 0.260610696402935d0 0.312347077040003d0 0.33023935501260d0
                                          0.312347077040003d0 0.260610696402935d0 0.180648160694857d0 0.081274388361574d0)))))))))

(defmacro gauss-legendre/n (f &key (n 1) (integrator 'gauss-legendre/4))
  "Create a lambda that integrates a region over `F`. The region will
be divided into `N` segments, and `INTEGRATOR` will be used to
integrate each segment. The lambda takes the lower and upper limits of
the region as parameters."
  `(let ((integrator (,integrator ,f)))
     (lambda (a b)
       (let ((h (/ (- b a) ,n)))
         (loop repeat ,n
            as l = a then (+ l h)
            as u = (+ a h) then (+ u h)
            summing (funcall integrator l u))))))

(defmacro gauss-legendre (f &key (integrator 'gauss-legendre/4) (tolerance +tolerance+) (limit +limit+))
  "Create a lambda that integrates a region over `F`.
`INTEGRATOR` specifies the algorithm to use. The region will be
subdivided into progressively smaller regions in order to converge to
the specified accuracy, within the specified iteration limit. The
lambda takes the lower and upper limits of the region as parameters."
  `(lambda (a b)
     (loop for x from 2 to (1- ,limit)
        as iny = (funcall (gauss-legendre/n ,f :n x :integrator ,integrator) a b) then inx
        as inx = (funcall (gauss-legendre/n ,f :n (1+ x) :integrator ,integrator) a b)
        as corr = (abs (/_ (- inx iny) inx))
        do (when (> ,tolerance corr) (return (values inx t x)))
        finally (return (values inx nil x)))))

(defmacro laguerre-quadrature/n (f n)
  "Create a lambda that integrates the region from 0 to infinity over `F`.
`N` specifies the order of the Laguerre polynomial to use for the
integration rule, which can be 2,4,6, or 8."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (flet ((lq (ti ci)
                (lambda ()
                  (loop for v in ti
                     for z in ci
                     summing (* z (funcall ,g v))))))
         (apply #'lq (cond ((eql ,n 2) '((.5857864376269050d0 3.414213562373095d0)
                                         (.8535533905932738d0 .1464466094067262d0)))
                           ((eql ,n 4) '((.3225476896193923d0 1.745761101158347d0
                                          4.536620296921128d0 9.395070912301133d0)
                                         (.6031541043416336d0 .3574186924377997d0
                                          .03888790851500538d0 .0005392947055613275d0)))
                           ((eql ,n 6) '((.2228466041792607d0 1.188932101672623d0
                                          2.992736326059314d0 5.775143569104511d0
                                          9.837467418382590d0 15.98287398060170d0)
                                         (.4589646739499636d0 .4170008307721210d0
                                          .1133733820740450d0 .01039919745314907d0
                                          .0002610172028149321d0 .0000008985479064296212d0)))
                           ((eql ,n 8) '((.1702796323051010d0 .9037017767993799d0
                                          2.251086629866131d0 4.266700170287659d0
                                          7.045905402393466d0 10.75851601018100d0
                                          15.74067864127800d0 22.86313173688926d0)
                                         (.3691885893416375d0 .4187867808143430d0
                                          .1757949866371718d0 .03334349226121565d0
                                          .002794536235225673d0 .00009076508773358213d0
                                          .0000008485746716272532d0 .000000001048001174871510d0)))))))))

(eval-when (:execute :load-toplevel :compile-toplevel)
  (defun defgaussian (n)
    `(defmacro ,(intern (format nil "GAUSS-LEGENDRE/~A" n)) (f)
       ,(format nil "Create a lambda that integrates a region over `F` using a Legendre polynomial~%of order ~A. The lambda takes the lower and upper limits of the region as parameters." n)
       `(gaussian-quadrature/n ,f ,,n)))
  (defun deflaguerre (n)
    `(defmacro ,(intern (format nil "GAUSS-LAGUERRE/~A" n)) (f)
       ,(format nil "Create a lambda that integrates the region from 0 to infinity over `F` using a Laguerre polynomial of order ~A." n)
       `(laguerre-quadrature/n ,f ,,n)))
  (loop for i from 2 to 9 do (eval (defgaussian i)))
  (loop for i from 2 to 8 by 2 do (eval (deflaguerre i))))

(defun running-std-dev ()
  "Create a lambda that can be called with successive values. The lambda returns the
standard deviation of the values it has been passed."
  (let ((sum1 0d0)
        (sum2 0d0)
        (n 0))
    (lambda (xi)
      (setq sum1 (+ sum1 (* xi xi)))
      (setq sum2 (+ sum2 xi))
      (incf n)
      (if (eq 1 n)
          0d0
          (sqrt (/ (- (* (/ 1d0 n) sum1) (expt (* (/ 1d0 n) sum2) 2))
                   (1- n)))))))

(defun standard-deviation (lx)
  "Compute the standard deviation of the values in the list `LX`."
  (let ((rsd (running-std-dev)))
    (loop for xi in lx
       as v = (funcall rsd xi)
       finally (return v))))

(defmacro running-monte-carlo (f a b)
  "Create a lambda that can be called with successive values ranging
randomly over the range from `A` to `B`. The lambda returns the
estimated integral of `F` between `A` and `B`, and the standard
deviation of the estimate."
  (let ((g (gensym)))
    `(let ((,g (genlambda/1 ,f)))
       (let ((n 0)
             (sd (running-std-dev))
             (c (- ,b ,a))
             (sum 0d0))
         (lambda (x)
           (let ((xi (funcall ,g x)))
             (incf n)
             (setq sum (+ sum xi))
             (values (* c (/ 1 n) sum)
                     (funcall sd xi))))))))

(defmacro monte-carlo (f a b &key (n 1000) (random (lambda () (random 1d0 (make-random-state t)))))
  "Use monte-carlo integration to estimate the integral of `F` from
`A` to `B`. Use `N` samples for the estimate. `RANDOM` should be a
lambda that returns a value between 0 and 1 every time it is called."
  `(flet ((xi ()
            (+ ,a (* (- ,b ,a) (funcall ,random)))))
     (let ((rmc (running-monte-carlo ,f ,a ,b)))
       (loop repeat (1- ,n)
          do (funcall rmc (xi)))
       (funcall rmc (xi))))) ; pass multiple value return from final rmc call
