{-# LANGUAGE BangPatterns #-}
module Main where
import Bhut (sun, earth, doUpdate, demo, mass, Body(..), position, velocity)
import Control.Concurrent
import Data.IORef
import Debug.Trace
import SDL hiding (trace)
import SDL.Vect(V2(..), V4(..))
import SDL.Primitive(fillCircle)
import Foreign.C.Types

import qualified Data.Text as T
import Control.Monad (unless, forever, void)

test_bodies = [ Body { mass = 10,
                       position = V2 100 100,
                       velocity = V2 0 0 }
              , Body { mass = 15,
                       position = V2 100 10,
                       velocity = V2 0 0 }
              , Body { mass = 10,
                       position = V2 500 500,
                       velocity = V2 0 0 }
              , Body { mass = 10,
                       position = V2 700 700,
                       velocity = V2 0 0 } ]

toCInt :: Rational -> CInt
toCInt r = (round . fromRational) r

adjustToOrigin :: Renderer -> V2 Rational -> IO (V2 Rational)
adjustToOrigin renderer (V2 x y) = do
  viewport <- get $ rendererLogicalSize renderer
  case viewport of
    Nothing -> pure $ (V2 x y)
    Just (V2 w h) ->
      do
        let !blah = trace "viewport" (V2 w h)
        return $ V2 (x + (toRational w)) (y + (toRational h))
  
renderBodies :: Renderer -> [Body] -> IO ()
renderBodies renderer bodies = do
  rendererDrawColor renderer $= V4 0 0 0 255 -- black
  clear renderer

  rendererDrawColor renderer $= V4 0 255 0 255 -- green
  let points = map (P . fmap toCInt . position) bodies
  putStrLn $ "points: " ++ show points  
  adjusted <-  pure -- $ mapM id
  -- (adjustToOrigin renderer)
                    $ map position bodies
  putStrLn $ "adjusted: " ++ (show (map (P . fmap toCInt) adjusted))
  mapM_ (\pos -> fillCircle renderer pos 4 (V4 0 255 0 255)) $ map (fmap toCInt) adjusted
  present renderer
  
main = do
  initializeAll
  window <- createWindow (T.pack "physics go brr") defaultWindow
  renderer <- createRenderer window (-1) defaultRenderer

  bodies <- newIORef test_bodies
  
  let loop = (do
                 pollEvents
                 bs <- readIORef bodies
                 putStrLn "next iter"
                 -- putStrLn $ show bs
                 renderBodies renderer bs
                 modifyIORef' bodies (doUpdate 20000)
                 loop)
  loop
  
  destroyWindow window
