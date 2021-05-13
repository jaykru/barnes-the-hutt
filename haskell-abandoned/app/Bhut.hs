{-# LANGUAGE NamedFieldPuns #-}
module Bhut (demo, earth, sun, Body(..), mass, position, velocity, doUpdate) where
import Data.Maybe
import SDL.Vect(V2(..))
  
type Vector = V2 Rational
data Body = Body
            {
              mass :: Rational
            , position :: Vector
            , velocity :: Vector
            } deriving (Show, Eq)

computeCenter :: [Body] -> Vector
computeCenter bodies =
  let totalMass = foldl1 (+) $ map (\body -> mass body) bodies
      weightedsum = foldl1 (+) $ map (\body -> (position body) * fromRational (mass body)) bodies
  in weightedsum / fromRational totalMass

type Extent = Maybe (Rational, Rational, Rational, Rational)

inExtentDec :: Body -> Extent -> Bool
inExtentDec body Nothing = False
inExtentDec Body{position= V2 x y} (Just(xmin, xmax, ymin, ymax)) = xmin <= x && x <= xmax && ymin <= y && y <= ymax

computeExtent :: [Body] -> Extent
computeExtent [] = Nothing
computeExtent (body:bodies) =
  let V2 x y = position body in
  case computeExtent bodies of
    Just (xmin,xmax,ymin,ymax) -> Just (min xmin x, max xmax x, min ymin y, max ymax y)
    Nothing -> Just (x, x, y, y)

data Quadtree = Quadtree
                {
                  body :: Maybe Body -- body stored here in the case of a leaf node
                , extent :: Extent
                , treemass :: Rational -- total mass of the tree in kilograms
                , treecenter :: Maybe Vector
                , q1 :: Maybe Quadtree
                , q2 :: Maybe Quadtree
                , q3 :: Maybe Quadtree
                , q4 :: Maybe Quadtree
                } deriving Show

emptyQuadtree = Quadtree { body = Nothing,
                           extent = Nothing,
                           treemass = 0,
                           treecenter = Nothing,
                           q1 = Nothing,
                           q2 = Nothing,
                           q3 = Nothing,
                           q4 = Nothing }

singletonQuadtree body = Quadtree { body = Just body,
                                    extent = computeExtent [body],
                                    treemass = mass body,
                                    treecenter = Just $ position body,
                                    q1 = Nothing,
                                    q2 = Nothing,
                                    q3 = Nothing,
                                    q4 = Nothing }


  
allBodies :: Quadtree -> [Body]
allBodies Quadtree { body, q1, q2, q3, q4 } =
  let rest = concat $ map allBodies $ catMaybes [q1, q2, q3, q4] in
    case body of Nothing -> []
                 Just body -> [body] ++ rest
  
buildQuadtree :: [Body] -> Quadtree
buildQuadtree [] = emptyQuadtree
buildQuadtree [body] = singletonQuadtree body
buildQuadtree bodies =
  let Just (xmin, xmax, ymin, ymax) = computeExtent bodies in
  let xavg = (xmin + xmax) / 2
      yavg = (ymin + ymax) / 2

      q1Extent = Just (xavg, xmax, yavg, ymax)
      q2Extent = Just (xmin, xavg, yavg, ymax)
      q3Extent = Just (xmin, xavg, ymin, yavg)
      q4Extent = Just (xavg, xmax, ymin, yavg)

      q1Bodies = [body | body <- bodies, body `inExtentDec` q1Extent]
      q2Bodies = [body | body <- bodies, body `inExtentDec` q2Extent]
      q3Bodies = [body | body <- bodies, body `inExtentDec` q3Extent]
      q4Bodies = [body | body <- bodies, body `inExtentDec` q4Extent]
  in
    Quadtree { body = Nothing,
               extent = Just (xmin, xmax, ymin, ymax),
               treemass = foldl1 (+) $ map mass bodies,
               treecenter = Just $ computeCenter bodies,
               q1 = Just $ buildQuadtree q1Bodies,
               q2 = Just $ buildQuadtree q2Bodies,
               q3 = Just $ buildQuadtree q3Bodies,
               q4 = Just $ buildQuadtree q4Bodies
               }

withoutBody :: Quadtree -> Body -> Quadtree
withoutBody tree body =
  let bodies = allBodies tree
      (ok, just_one) = break (== body) bodies
  in buildQuadtree (ok ++ (drop 1 just_one))

without :: Eq a => a -> [a] -> [a]
without a as =
  let (ok, just_one) = break (== a) as
  in ok ++ drop 1 just_one

-- fixed accuracy parameter, roughly 1
theta :: Rational
theta = 1

euclideanDist :: Vector -> Vector -> Rational
euclideanDist (V2 x1 y1) (V2 x2 y2) =
  toRational $ sqrt $ fromRational $ (x1 - x2)^2 + (y1 - y2)^2

magnitude :: Vector -> Rational
magnitude = euclideanDist $ V2 0 0

-- computes force of gravity between two bodies as a vector indicating
-- the force applied on object 1 by object 2
newtonianForce :: Body -> Body -> Vector
newtonianForce Body{mass=m1, position=pos1} Body{mass=m2, position=pos2} =
  unit * fromRational ((gConst * m1*m2) / r^2) -- (N m^2 / kg^2) * kg^2 / m^2 = N
  where gConst = 6.674 * (toRational $ 10**(-11)) -- (N m^2) / kg^2
        r = euclideanDist pos1 pos2 -- m
        unit = (pos2 - pos1) * fromRational (magnitude $ pos2 - pos1)

computeForce :: Body -> Maybe Quadtree -> Vector
computeForce Body{mass, position, velocity} tree =
  case tree of
    Nothing -> V2 0 0
    Just tree ->
      case extent tree of
        Nothing -> V2 0 0
        Just (xmin, xmax, ymin, ymax) -> 
          let Just center = treecenter tree
              body = Body{mass=mass, position=position, velocity=velocity}
              diameter = euclideanDist center position
              -- TODO: is this actually the right calculation?
              cellLength = max (abs $ xmin - xmax) (abs $ ymin - ymax) 
          in
            if cellLength / diameter < theta then
              newtonianForce body Body{mass = treemass tree, position = center}
            else
              foldl1 (+) $ [ computeForce body subtree
                           | subtree <- [q1 tree, q2 tree, q3 tree, q4 tree]]

nextVelocity :: Body -> Vector -> Rational -> Vector
nextVelocity body force time =
  let v0 = velocity body
      accel = force / fromRational (mass body)
  in
    v0 + accel * fromRational time

rat :: Rational -> V2 Rational
rat = fromRational

nextPosition :: Body -> Vector -> Rational -> Vector
nextPosition body force time =
  let pos0 = position body
      accel = force / rat (mass body)
      v0 = velocity body
  in
    pos0 + v0* rat time + accel * rat (1/2 * time^2)

doUpdate :: Rational -> [Body] -> [Body]
doUpdate time bodies =
  let aux (b:bs) = 
        let force = computeForce b (Just $ buildQuadtree (without b bodies))
            nextVel = nextVelocity b force time
            nextPos = nextPosition b force time
            nextBody = b { position = nextPos, velocity = nextVel }
        in
          nextBody:(aux bs)
      aux [] = []
  in aux bodies

earth = Body { mass = 5.972 * 10^24
             ,
               position = V2 (toRational 1.49*10^11) (toRational 1.49*10^11),
               velocity = V2 0 0 }
  
sun = Body { mass = 2 * 10^30
           ,
             position = V2 0 0,
             velocity = V2 0 0 }

demo =
  let qt = buildQuadtree [sun] -- NOTE: we probably need to add a
                               -- function to remove the current node
                               -- from consideration for use with the
                               -- force computation

  in
    do putStrLn $ show $ buildQuadtree [sun,earth]
       putStrLn $ "Total force on the earth is " ++ (show $ computeForce earth (Just qt)) ++ "N"
       putStrLn $ "The newtonian calculation yields " ++ (show $ newtonianForce earth sun) ++ "N"
       putStrLn $ show $ doUpdate 1 [earth,sun]
  
