(in-package :cl-user)

(defpackage :numerical.roots
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+)
  (:export :bisect
           :bisect+
           :newton-raphson
           :newton-raphson+
           :nr/bisect
           :nr/bisect+
           :root-brackets
           :root-brackets+
           :false-position
           :false-position+
           :secant/bisect
           :secant/bisect+
           :quadratic/bisect
           :quadratic/bisect+))

(in-package :numerical.roots)

(defun bisect (f &key (tolerance +tolerance+) (limit +limit+))
  (labels ((rec (left right &optional (cnt 0) (fleft (funcall f left)) (fright (funcall f right)))
             (let* ((middle (* 0.5d0 (+ left right)))
                    (fmid (funcall f middle)))
               (if (or
                    (> (* fleft fright) 0d0) ; not bracketed
                    (>= cnt limit))
                   (values middle nil cnt)
                   (if (< (abs (/_ (- right left) middle)) tolerance)
                       (values middle t cnt)
                       (if (< (* fleft fmid) 0)
                           (rec left middle (1+ cnt) fleft fmid)
                           (rec middle right (1+ cnt) fmid fright)))))))
    #'rec))

(defmacro bisect+ ((x &key (tolerance +tolerance+) (limit +limit+)) &rest rest)
  `(bisect (lambda (,x) ,@rest) :tolerance ,tolerance :limit ,limit))

(defun newton-raphson (f df &key (tolerance +tolerance+) (limit +limit+))
  (labels ((rec (guess &optional (cnt 0))
             (if (>= cnt limit)
                 (values guess nil limit)
                 (let ((delta (/_ (funcall f guess) (funcall df guess) -1d0)))
                   (if (<= (abs (/_ delta (+ delta guess))) tolerance)
                       (values (+ guess delta) t cnt)
                       (rec (+ guess delta) (1+ cnt)))))))
    #'rec))

(defmacro newton-raphson+ ((x &key (tolerance +tolerance+) (limit +limit+)) f df)
  `(newton-raphson (lambda (,x) ,f) (lambda (,x) ,df) :tolerance ,tolerance :limit ,limit))

(defun nr/bisect (f df &key (tolerance +tolerance+) (limit +limit+))
  (labels ((rec (left right &optional
                      (cnt 0)
                      (fleft (funcall f left))
                      (fright (funcall f right))
                      (best (if (<= (abs fleft) (abs fright)) left right))
                      (fbest (funcall f best))
                      (dfbest (funcall df best))
                      (nr (lambda () (let ((delta (/_ fbest dfbest -1d0)))
                                       (values delta (+ best delta)))))
                      (bisect (lambda () (values (/_ (- right left) 2d0)
                                                 (/_ (+ right left) 2d0)))))
             (if (or
                  (> (* fleft fright) 0d0) ; not bracketed
                  (> cnt limit)) ; fail to converge
                 (values best nil cnt)
                 (multiple-value-bind (delta best) (funcall (if (<= (*
                                                                     (- (* dfbest (- best left)) fbest)
                                                                     (- (* dfbest (- best right)) fbest))
                                                                    0d0)
                                                                nr bisect))
                   (if (< (abs (/_ delta best)) tolerance)
                       (values best t cnt)
                       (if (>= (* (funcall f best) fleft) 0d0)
                           (rec best right (1+ cnt))
                           (rec left best (1+ cnt))))))))
    #'rec))

(defmacro nr/bisect+ ((x &key (tolerance +tolerance+) (limit +limit+)) f df)
  `(nr/bisect (lambda (,x) ,f) (lambda (,x) ,df) :tolerance ,tolerance :limit ,limit))

(defun root-brackets (f &key (steps 100))
  (lambda (left right)
    (loop
       for px = left then x
       for x = left then (+ x (/ (- right left) steps))
       when (not (equal (plusp (funcall f x))
                        (plusp (funcall f px))))
       collect (list (coerce px 'double-float) (coerce x 'double-float)) into res
       until (>= x right)
       finally (return res))))

(defmacro root-brackets+ ((x &key (steps 100)) &rest rest)
  `(root-brackets (lambda (,x) ,@rest) :steps ,steps))

(defun false-position (f &key (tolerance +tolerance+) (limit +limit+))
  (labels ((rec (left right &optional
                      (cnt 0)
                      (prev right)
                      (fleft (funcall f left))
                      (fright (funcall f right)))
             (if (> cnt limit)
                 (values prev nil cnt)
                 (let* ((guess (/_ (- (* left fright) (* 1d0 right fleft))
                                   (- fright fleft)))
                        (fguess (funcall f guess))
                        (err (abs (/_ (- guess prev) guess))))
                   (if (<= err tolerance)
                       (values guess t cnt)
                       (if (< (* fleft fguess) 0d0)
                           (rec left guess (1+ cnt) guess fleft fguess)
                           (rec guess right (1+ cnt) guess fguess fright)))))))
    #'rec))

(defmacro false-position+ ((x &key (tolerance +tolerance+) (limit +limit+)) &rest rest)
  `(false-position (lambda (,x) ,@rest) :tolerance ,tolerance :limit ,limit))

(defun secant/bisect (f &key (tolerance +tolerance+) (limit +limit+))
  (let ((bistepper (bisect f :tolerance tolerance :limit 1)))
    (labels ((secant (xi-1 xi)
               (let ((fxi-1 (funcall f xi-1))
                     (fxi (funcall f xi)))
                 (- xi (* 1d0 fxi (/ (- xi xi-1)
                                     (- fxi fxi-1))))))
             (rec (left right &optional
                        (cnt 0)
                        (fleft (funcall f left)) 
                        (fright (funcall f right))
                        (best (if (<= (abs fleft) (abs fright)) left right))
                        (xi-1 (funcall bistepper left right 1))
                        (xi (funcall bistepper left right 0)))
               (if (> (* fleft fright) 0d0)
                   (values best nil cnt)
                   (if (> cnt limit)
                       (funcall (bisect f :tolerance tolerance :limit limit) left right)
                       (let* ((guess (secant xi-1 xi))
                              (fguess (funcall f guess))
                              (err (abs (/_ (- guess xi) guess))))
                         (if (<= err tolerance)
                             (values guess t cnt)
                             (if (<= (* fleft fguess) 0d0)
                                 (rec left guess (1+ cnt) fleft fguess guess xi guess)
                                 (rec guess right (1+ cnt) fguess fright guess xi guess))))))))
      #'rec)))

(defmacro secant/bisect+ ((x &key (tolerance +tolerance+) (limit +limit+)) &rest rest)
  `(secant/bisect (lambda (,x) ,@rest) :tolerance ,tolerance :limit ,limit))

(defun quadratic/bisect (f &key (tolerance +tolerance+) (limit +limit+))
  (let ((bistepper (bisect f :tolerance tolerance :limit 2)))
    (labels ((quadratic (x0 x1 x2)
               (if (or
                    (equal x0 x1)
                    (equal x1 x2)
                    (equal x0 x2))
                   x2
                   (let* ((fx0 (funcall f x0))
                          (fx1 (funcall f x1))
                          (fx2 (funcall f x2))
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
                        (fleft (funcall f left)) 
                        (fright (funcall f right))
                        (best (if (<= (abs fleft) (abs fright)) left right))
                        (xi-2 (funcall bistepper left right 2))
                        (xi-1 (funcall bistepper left right 1))
                        (xi (funcall bistepper left right 0)))
               (if (> (* fleft fright) 0d0)
                   (values best nil cnt)
                   (if (> cnt limit)
                       (funcall (bisect f :tolerance tolerance :limit limit) left right)
                       (let* ((guess (quadratic xi-2 xi-1 xi))
                              (fguess (funcall f guess))
                              (err (abs (/_ (- guess xi) guess))))
                         (if (<= err tolerance)
                             (values guess t cnt)
                             (if (<= (* fleft fguess) 0d0)
                                 (rec left guess (1+ cnt) fleft fguess guess xi-1 xi guess)
                                 (rec guess right (1+ cnt) fguess fright guess xi-1 xi guess))))))))
      #'rec)))

(defmacro quadratic/bisect+ ((x &key (tolerance +tolerance+) (limit +limit+)) &rest rest)
  `(quadratic/bisect (lambda (,x) ,@rest) :tolerance ,tolerance :limit ,limit))
