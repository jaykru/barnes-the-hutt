(declaim (optimize (speed 0) (space 0) (debug 3)))
(ql:quickload '(:alexandria
                :array-operations
                :clack
                :woo
                :cl-json
                :trivial-open-browser
                :websocket-driver
                :ieee-floats
                :lparallel
                :cl-cpus))
(ignore-errors (ql-dist:install-dist "http://dist.ultralisp.org/"
                                     :prompt nil))
(ql:quickload :queue)
(use-package :queue)

(defclass body ()
  ((mass :initarg :mass)                ; in SI kilograms
   (position :initarg :position)        ; a vector #(x y)
   (velocity :initarg :velocity)))      ; in m/s


(defparameter *test-bodies*
  (loop for x from 0 to 150
        append (loop for y from 0 to 100
                     collect (make-instance 'body :position `#(,(coerce x 'float) ,(coerce y 'float))
                                                  :mass 1
                                                  :velocity #(0 0)))))
(defmethod print-object ((b body) stream)
  (with-slots (mass position velocity) b
    (let ((x (elt position 0))
          (y (elt position 1))
          (vx (elt velocity 0))
          (vy (elt velocity 1)))
      (print-unreadable-object (b stream :type t)
        (format stream "{ [~a,~a]m ; [~a,~a]m/s; ~akg }" x y vx vy mass)))))

;; decides whether a certain body falls within an extent
(defun inextp (body extent)
  "Decides wheter a body falls within an extent"
  (with-slots (position) body
    (let ((xmin (car extent))
          (xmax (cadr extent))
          (ymin (caddr extent))
          (ymax (cadddr extent))

          (x (elt position 0))
          (y (elt position 1)))
      (and (and (<= xmin x) (<= x xmax))
           (and (<= ymin y) (<= y ymax))))))

(defun v+ (&rest vecs)
  "Adds two vectors"
  (if (not vecs)
      #(0 0)
      (let* ((v1 (car vecs))
             (v2 (cadr vecs))
             (v1x (elt v1 0))
             (v1y (elt v1 1))
             (v2x (elt v2 0))
             (v2y (elt v2 1)))
        `#(,(+ v1x v2x) ,(+ v1y v2y)))))

(defun v/ (v1 k)
  "Divides a vector by a non-zero scalar"
  (map 'vector #'(lambda (x) (/ x k)) v1))

(defun v* (v1 k)
  "Scales a vector by a scalar"
  (map 'vector #'(lambda (x) (* x k)) v1))

(defun v- (v1 v2)
  "Subtracts two vectors"
  (v+ v1 (v* v2 -1)))

(defun compcenter (bodies)
  "Computes the center of mass of some bodies"
  (if (not bodies)
      (error "tried to compute center of no bodies")
      (let ((M (compmass bodies)))
        (v/ (reduce #'v+
                    (mapcar #'(lambda (body)
                                (with-slots (mass position) body
                                  (v* position mass)))
                            bodies))
            M))))

(defun compextent (bodies)
  "Compute the geometric extent of some bodies"
  (cond
    ((equal bodies nil)
     nil)
    ('t
     (loop for body in bodies
           maximizing (elt (slot-value body 'position) 0) into xmax
           minimizing (elt (slot-value body 'position) 0) into xmin
           maximizing (elt (slot-value body 'position) 1) into ymax
           minimizing (elt (slot-value body 'position) 1) into ymin
           finally (return `(,xmin ,xmax ,ymin ,ymax))))))

(defclass qtree ()
  ((body :initarg :body
         :initform nil) ;; the body contained in this node, nil unless the
   ;; quadtree is a leaf

   (extent :initarg :extent
           :initform nil)             ; a list '(xmin xmax ymin max) which keeps
                                        ; track of the area a given qtree node
                                        ; covers

   ;; the mass of all the objects contained within the extent of this tree
   (treemass
    :initarg :treemass
    :initform nil)
   ;; in the case of a (non-empty) leaf, this is just the object mass

   ;; center of mass of the extent represented by this tree
   (treecenter
    :initarg :treecenter
    :initform nil)

   ;; quadrants
   (q1 :initarg :q1 :initform nil)
   (q2 :initarg :q2 :initform nil)
   (q3 :initarg :q3 :initform nil)
   (q4 :initarg :q4 :initform nil)))

;; makes a qtree with an initial extent
(defun make-qtree (start-ext)
  (make-instance 'qtree :extent start-ext))

(defun leafp (tree)
  (equal (with-slots (q1 q2 q3 q4) tree
           `(,q1 ,q2 ,q3 ,q4))
         `(nil nil nil nil)))

(defun emptyp (tree)
  (or (not tree)
      (and (not (slot-value tree 'body))
           (leafp tree))))


(defmacro insertff (body trees)
  "Insert a body into the first-fitting subtree. This is a macro so that SBCL can
    optimize the tail-recursion."
  `(loop for tree in ,trees do
    (with-slots (extent) tree
      (if (inextp ,body extent)
          (return (insert-body tree ,body))))))

(defun allbodies (tree)
  (if (not tree)
      nil
      (with-slots (body q1 q2 q3 q4) tree
        (remove-if #'(lambda (b) (equal nil b)) (cons body (append (allbodies q1)
                                                              (allbodies q2)
                                                              (allbodies q3)
                                                              (allbodies q4)))))))

(defun quarter-extent (extent)
  (destructuring-bind (xmin xmax ymin ymax) extent
    (let* ((xavg (coerce  (/ (+ xmin xmax) 2) 'double-float))
           (yavg (coerce (/ (+ ymin ymax) 2) 'double-float))

           (q1ext `(,xavg ,xmax ,yavg ,ymax))
           (q2ext `(,xmin ,xavg ,yavg ,ymax))
           (q3ext `(,xmin ,xavg ,ymin ,yavg))
           (q4ext `(,xavg ,xmax ,ymin ,yavg)))
      `(,q1ext ,q2ext ,q3ext ,q4ext))))

(defmethod quarter ((tree qtree))
  "Splits a leaf into quadrants, placing its body (if there is one) into the
   appropriate new quadrant."
  (with-slots (extent q1 q2 q3 q4 body) tree
    (if (and (leafp tree)               ; leaf
             extent)                   ; extentful
        (destructuring-bind (q1ext q2ext q3ext q4ext) (quarter-extent extent)
          (let* ((q1p (make-qtree q1ext))
                 (q2p (make-qtree q2ext))
                 (q3p (make-qtree q3ext))
                 (q4p (make-qtree q4ext)))
            (progn (setf q1 q1p)
                   (setf q2 q2p)
                   (setf q3 q3p)
                   (setf q4 q4p)
                   (if body
                     (let ((oldbody body))
                       (setf body nil)
                       (insertff oldbody `(,q1 ,q2 ,q3 ,q4)))))))
        (error "We only quarter leaves with an extent"))))

(defun compmass (bodies)
  (reduce #'+ (mapcar (lambda (b) (slot-value b 'mass)) bodies)))

(defmethod insert-body ((tree qtree) b)
  (with-slots (body extent treemass treecenter q1 q2 q3 q4) tree
    (if (inextp b extent)
        (if (not body)
            (if (leafp tree)            ; just fill the empty leaf
                (progn
                  (setf body b)
                  (setf treemass (slot-value b 'mass))
                  (setf treecenter (slot-value b 'position)))
                                        ; here we have a node already split into
                                        ; quadtrees, so we just need to insert b
                                        ; in to the first fitting subtree and
                                        ; update all the tree metrics
                                        ; appropriately
                (progn (setf treecenter (compcenter (cons b (allbodies tree))))
                       (setf treemass (compmass (cons b (allbodies tree))))
                       (insertff b `(,q1 ,q2 ,q3 ,q4))))
                                        ; full leaf: in this case we need to
                                        ; split the leaf, as we only allow one
                                        ; body per leaf
            (if (not (equalp body b))
                (progn
                  (quarter tree)
                  (setf treecenter (compcenter (cons b (allbodies tree))))
                  (setf treemass (compmass (cons b (allbodies tree))))
                  (insertff b `(,q1 ,q2 ,q3 ,q4)))))
        (error "Tried inserting a body into a tree with wrong extent"))))

(defun eq-body (b1 b2)
  (if (and (not b1) (not b2))
      't
      (and (equalp (slot-value b1 'position) (slot-value b2 'position))
           (equalp (slot-value b1 'mass) (slot-value b2 'mass)))))

(defun build-qtree (bodies)
  (let ((init-tree (make-qtree (compextent bodies))))
    (progn
      (loop for bod in bodies do
        (insert-body init-tree bod))
      init-tree)))

(defun euclidean-dist (v1 v2)
  (let ((x1 (elt v1 0))
        (x2 (elt v2 0))
        (y1 (elt v1 1))
        (y2 (elt v2 1)))
    (sqrt (+ (expt (- x1 x2) 2) (expt (- y1 y2) 2)))))

(defun magnitude (v)
  (euclidean-dist v #(0 0)))

(defun newtonian-force (b1 b2)
  (if (or (equalp (slot-value b1 'position) (slot-value b2 'position)) (eq-body b1 b2))
      #(0 0)
      (let* ((p1 (slot-value b1 'position))
             (p2 (slot-value b2 'position))
             (m1 (slot-value b1 'mass))
             (m2 (slot-value b2 'mass))
             (unit (v* (v- p1 p2) (- (magnitude (v- p2 p1)))))
             (r (euclidean-dist p1 p2))
             (gconst (* 6.674 (expt 10 -11))))
        (v* unit (/ (* gconst m1 m2) (expt r 2))))))

(defparameter *theta* 1)

(defmethod centroid ((tree qtree))
  (with-slots (treecenter treemass) tree
    (make-instance 'body :mass treemass :position treecenter :velocity #(0 0))))

(defmethod subtrees ((tree qtree))
  (with-slots (q1 q2 q3 q4) tree
    `(,q1 ,q2 ,q3 ,q4)))

(defun comp-force (body tree)
  (with-slots (extent treecenter treemass q1 q2 q3 q4) tree
    (if (or (emptyp tree) (not extent))
        #(0 0)
        (if (not treecenter)
            (error "Trying to compute force due to tree with no center")
            (let* ((xmin (car extent))
                   (xmax (cadr extent))
                   (ymin (caddr extent))
                   (ymax (cadddr extent))
                   (diameter (euclidean-dist treecenter (slot-value body 'position)))
                   (cell-length (max (abs (- xmin xmax)) (abs (- ymin ymax)))))
              (if (or (equalp diameter 0.0)
                      (< (/ cell-length diameter) *theta*))
                  (newtonian-force body (centroid tree))
                  (reduce #'v+ (mapcar #'(lambda (subtree) (comp-force body subtree))
                                       (subtrees tree)))))))))

(defun next-velocity (body force time)
  (let* ((v0 (slot-value body 'velocity))
         (m (slot-value body 'mass))
         (accel (v/ force m)))
    (v+ v0 (v* accel time))))

(defun next-position (body force time)
  (let ((pos0 (slot-value body 'position))
        (accel (v/ force (slot-value body 'mass)))
        (v0 (slot-value body 'velocity)))
    (v+ pos0
        (v+ (v* v0 time) (v* accel (* 1/2 (expt time 2)))))))

(defun update-bodies (time bodies)
  (let ((tree (par-build-qtree bodies)))
    (loop for body in bodies do
      (let ((force (comp-force body tree)))
        (with-slots (velocity position) body
          (progn
            (setf position (next-position body force time))
            (setf velocity (next-velocity body force time))))))))

;; parallel stuff
(defun bits-of (n)
  "Computes the bits ***in LSB-first order*** of n in O(log(n)) time"
  (let* ((len (if (equalp n 0)
                  1
                  (1+ (floor (log n 2)))))
         (a (make-array len :initial-element 0)))
    (loop for i from 0 to (1- len) do
      (setf (aref a i)
            (if (equalp 0 (logand n (ash 1 i)))
                0
                1)))
    a))

(defun counting-sort (seq range &key (metric #'identity))
  "Sorts arr whose elements are in [0,range]. Optionally uses the function `metric`
   to determine values for ordering the value of range has different semantics
   when used with non-identity key"
  (let ((counts (make-array (1+ range) :initial-element 0)))
    (loop for j from 0 below (length seq) do
      (let ((key (funcall metric (nth j seq))))
        (setf (aref counts key) (1+ (aref counts key)))))
    ;; counts[(key v)] now contains the number of elements e in seq with (key e) equal to (key v)
    (loop for i from 1 to range do
      (setf (aref counts i) (+ (aref counts i) (aref counts (1- i)))))
    ;; counts[(metric v)] now contains the number of elements in seq with metric leq to (metric v)
    (let ((out (make-array (length seq) :initial-element 0)))
      (loop for j from (1- (length seq)) downto 0 do
        (setf (aref out (1- (aref counts (funcall metric (nth j seq))))) (nth j seq))
        (setf (aref counts (funcall metric (nth j seq))) (1- (aref counts (funcall metric (nth j seq))))))
      out)))

(defun left-pad-nth-bit (bits n)
  "Takes an array bits and gives the n-th bit, from right to left; giving 0 if
   we have exceeded the length of the bits."
  (handler-case (aref bits n)
    (error (c)
      (declare (ignore c))
      0)))

(defun radix-sort (seq &key (metric #'identity))
  "Sorts seq a list a natural numbers by using counting sort on each bit"
  (let* ((N (loop for v in seq
                  maximize (if (equalp (funcall metric v) 0)
                               1
                               (1+ (floor (log (funcall metric v) 2)))) into num_bits
                  finally (return num_bits)))
         (as-bits (make-hash-table))
         (new-seq (loop for v in seq
                        collect `(,v . ,(funcall metric v)))))
    (loop for key in (mapcar metric seq)
          do (setf (gethash key as-bits) (bits-of key)))
    (loop for d from 0 to N
          do (setf new-seq (map 'list #'identity
                                (counting-sort
                                 new-seq
                                 1
                                 :metric (lambda (b)
                                           (left-pad-nth-bit (gethash (cdr b) as-bits) d))))))
    (mapcar #'car new-seq)))

(defun intersperse-hack (x y)
  "intersperses the bits of x and y by some fancy bit-bashing

   borrowed from http://www.graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
   return by date: Common Lisp EOL"
  (let ((B '(#x55555555 #x33333333 #x0F0F0F0F #x00FF00FF))
        (S '(1 2 4 8))
        (x x)
        (y y))
    (setf x (logand (logior x (ash x (nth 3 S)))
                    (nth 3 B)))
    (setf x (logand (logior x (ash x (nth 2 S)))
                    (nth 2 B)))
    (setf x (logand (logior x (ash x (nth 1 S)))
                    (nth 1 B)))
    (setf x (logand (logior x (ash x (nth 0 S)))
                    (nth 0 B)))

    (setf y (logand (logior y (ash y (nth 3 S)))
                    (nth 3 B)))
    (setf y (logand (logior y (ash y (nth 2 S)))
                    (nth 2 B)))
    (setf y (logand (logior y (ash y (nth 1 S)))
                    (nth 1 B)))
    (setf y (logand (logior y (ash y (nth 0 S)))
                    (nth 0 B)))
    (logior x (ash y 1))))

(defun morton-enc (body)
  "Produces a body's Morton code giving a linear encoding of its position in R^2"
  (with-slots (position) body
    (let ((x (coerce (elt position 0) 'float))
          (y (coerce (elt position 1) 'float)))
      (intersperse-hack (ieee-floats:encode-float64 y) (ieee-floats:encode-float64 x)))))

(defun morton-sort (bodies)
  (radix-sort bodies :metric #'(lambda (x) (morton-enc x))))

(defun schedule (bodies num-workers)
  (let* ((q (make-queue))
         (top (make-qtree (compextent bodies))))
    ;; start by morton-sorting the bodies so we ensure we get a good quadtree
    (enqueue q `(,bodies ,top ,`(0 . ,(1- (length bodies)))))
    (loop until (or
                 (< (length (caar (slot-value q 'queue::head))) 4)
                 (> (+ 4 (queue-size q)) num-workers)) do
                   (let* ((first (dequeue q))
                          (bs (car first))
                          (start-idx (caaddr first))
                          (width (ceiling (/ (length bs) 4)))

                          (tree (cadr first))

                          ;; quarter the bodies for the head of the queue and make a subtree
                          ;; skeleton for each quadrant
                          (chunks (loop for i from 0 to 3
                                        collect
                                        (let ((subbs (subseq bs
                                                             (* i width)
                                                             (min (length bs) (* (1+ i) width)))))
                                          `(,subbs ,(make-qtree (compextent subbs))
                                                   ,`(,(+ start-idx (* i width))
                                                      . ,(+ start-idx (min (1- (length bs))
                                                                           (1- (* (1+ i) width))))))))))
                     (loop for chunk in chunks
                           do (enqueue q chunk))
                     (with-slots (q1 q2 q3 q4 treemass treecenter extent) tree
                       (setf q1 (cadr (nth 0 chunks)))
                       (setf q2 (cadr (nth 1 chunks)))
                       (setf q3 (cadr (nth 2 chunks)))
                       (setf q4 (cadr (nth 3 chunks)))
                       (setf treemass (compmass bs))
                       (setf treecenter (compcenter bs))
                       (setf extent (compextent bs)))))
    (cons top (queue:queue-to-list q))))

(defun finish-build (skeleton my-range bodies)
  (destructuring-bind (lower . upper) my-range
    (loop for i from lower to upper do
      (insert-body skeleton (aref bodies i)))))

(setf lparallel:*kernel* (lparallel:make-kernel (cl-cpus:get-number-of-processors)))
(defun par-build-qtree (bodies)
  (let ((sorted (morton-sort bodies)))
    (destructuring-bind (top . chunks) (schedule sorted (cl-cpus:get-number-of-processors))
      (let ((body-array (map 'vector #'identity sorted)))
        (lparallel:pmap nil
                        #'(lambda (chunk)
                            (destructuring-bind (bs tree range) chunk
                              (declare (ignore bs))
                              (finish-build tree range body-array)))
                        chunks))
      top)))

(format t "parallel bench")
(time (par-build-qtree *test-bodies*))

(format t "sequential bench")
(time (build-qtree *test-bodies*))

;;; webapp
(defvar *conbods* (make-hash-table))

(defun handle-new-connection (con)
  (progn (setf (gethash con *conbods*) nil)))

(defun handle-close-connection (connection)
  (remhash connection *conbods*))

(defparameter *speed* 10000)
(defun handle-command (con message)
  (let ((m (json:decode-json-from-string message)))
    (let ((type (cdr (assoc :type m)))
          (bodies (cdr (assoc :data m))))
      (alexandria:switch (type :test #'equalp)
                         ("poll" (let ((old_bodies (gethash con *conbods*)))
                                   (if old_bodies
                                       (progn
                                         (update-bodies *speed* old_bodies)
                                         (websocket-driver:send
                                          con
                                          (json:encode-json-to-string
                                           (make-msg "update" (gethash con *conbods*))))))))
                         ("stop"
                          (setf (gethash con *conbods*) nil)
                          (websocket-driver:send
                           con
                           (json:encode-json-to-string
                            (make-msg "stop-ok" (gethash con *conbods*)))))
                         ("start"
                          (if bodies
                              (let* ((pos-seen (make-hash-table))
                                     (bods nil))
                                (loop for body in bodies
                                      do
                                         (let* ((x (cadr (assoc :position body)))
                                                (y (caddr (assoc :position body)))
                                                (pos `#(,x ,y))
                                                (mass (cdr (assoc :mass body)))
                                                (v1 (cadr (assoc :velocity body)))
                                                (v2 (caddr (assoc :velocity body)))
                                                (vel `#(,v1 ,v2)))
                                           (if (not (gethash pos pos-seen))
                                               (progn
                                                 (setf bods (cons
                                                             (make-instance 'body :position pos :mass mass :velocity vel) bods))))))
                                (websocket-driver:send con (json:encode-json-to-string (make-msg "ack" nil)))
                                (setf (gethash con *conbods*) bods))
                              (setf (gethash con *conbods*) nil)))
                         (otherwise (error "malformed message received"))))))

(defclass msg ()
  ((type :initarg :type :initform "update") ; type of websocket message
   (data :initarg :data :initform nil)))

(defun make-msg (type data)
  (make-instance 'msg :type type :data data))

(defun bhut-server (env)
  (let ((ws (websocket-driver:make-server env)))
    (websocket-driver:on :open ws
                         (lambda () (handle-new-connection ws)))

    (websocket-driver:on :message ws
                         (lambda (msg) (handle-command ws msg)))

    (websocket-driver:on :close ws
                         (lambda (&key code reason)
                           (declare (ignore code reason))
                           (handle-close-connection ws)))
    (lambda (responder)
      (declare (ignore responder))
      (websocket-driver:start-connection ws))))

(defparameter *html*
  (with-open-file (html (merge-pathnames *default-pathname-defaults* "bhut.html"))
    (alexandria:read-file-into-string html)))

(defun client-server (env)
  (declare (ignore env))
  `(200 (:content-type "text/html")
        (,(with-open-file (html (merge-pathnames *default-pathname-defaults* "bhut.html"))
            (alexandria:read-file-into-string html)))))

(clack:clackup #'client-server :port 8080)
(clack:clackup #'bhut-server :server 'woo :port 5000)
(trivial-open-browser:open-browser "http://localhost:8080")
