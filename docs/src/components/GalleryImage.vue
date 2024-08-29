<!-- Adding the thumbnail in src fails locally makes it that the page has to be reloaded 2x. Instead we use links.  -->
<script setup lang="ts">
import { withBase } from 'vitepress'
export interface Props {
  href: string;
  src: string;
  caption: string;
  desc: string;
}

defineProps<Props>();
</script>

<template>
  <div class="img-box">
    <a :href="href">
      <img :src="withBase(src)" height="150px" alt="">
      <div class="transparent-box1">
        <div class="caption">
          <h3>{{ caption }}</h3>
        </div>
      </div>
      <div class="transparent-box2">
        <div class="subcaption">
          <p class="opacity-low">{{ desc }}</p>
        </div>
      </div>
    </a>
  </div>
</template>

<style scoped>
.img-box {
  box-sizing: content-box;
  border-radius: 14px;
  margin: 10px;
  height: 190px;
  width: 140px;
  overflow: hidden;
  display: inline-block;
  color: white;
  position: relative;
  background-color: transparent;
  border: 2px solid var(--vp-c-bg-alt);
}

.img-box img {
  height: 100%;
  width: 100%;
  object-fit: cover;
  opacity: 0.4;
  transition: transform 0.3s ease, opacity 0.3s ease;
}

.caption {
  position: absolute;
  bottom: 30px;
  color: var(--vp-c-text-1);
  left: 10px;
  opacity: 1;
  transition: transform 0.3s ease, opacity 0.3s ease;
}

.subcaption {
  position: absolute;
  bottom: 5px;
  color: var(--vp-c-text-1);
  left: 10px;
  opacity: 0;
  transition: transform 0.3s ease, opacity 0.3s ease;
}

.subcaption p {
  line-height: 1.2;
  text-align: left;
  font-size: 10pt;
}

.transparent-box1 {
  height: 125px;
  width: 135px;
  background-color: transparent;
  position: absolute;
  top: 0;
  left: 0;
}

.transparent-box2 {
  height: 60px;
  width: 135px;
  background-color: transparent;
  position: absolute;
  top: 160px;
  left: 0;
}

.img-box:hover img {
  transform: scale(1.1);
  opacity: 0;
}

.img-box:hover .caption {
  transform: translateY(-20px);
  opacity: 1;
}

.img-box:hover .subcaption {
  transform: translateY(-20px);
  opacity: 1;
}

.img-box:hover {
  border: 2px solid var(--vp-c-brand-light);
  cursor: pointer;
}

.opacity-low {
  opacity: 0.85;
}
</style>